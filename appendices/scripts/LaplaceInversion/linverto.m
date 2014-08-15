function out = linvert(data, K, tau, varargin)
% This function takes a data structure (x, y, z) and calculates its 1D
% Laplace inversion with Tikhonov inversion. This is an ill-posed problem,
% and boils down to executing the minimization:
%
% Paul J. Ganssle, U.C. Berkeley, Pines Lab, June 2012
% This function largely implements the algorithm detailed in the paper:
%
% Venkataramanan, L., Song, Y.-Q., Hurlimann, M, D.,
%	"Solving Fredholm Integrals of the First Kind with Tensor Product
%		Structure in 2 and 2.5 Dimensions",
%
%		IEEE Transactions on Signal Processing, Vol. 50, No. 5 May 2002
%
% Released free for use, modification and non-commercial distribution, with
% credit.
%
% Inputs:
% data -> Struct, can be created from make_data_struct()
%   data.x = First dimension (e.g. time in indirect dimension 1) [size nx1]
%   data.y = Measurement response [size nx1]
%
% K    -> Kernel function. Can be anonymous function or function handles
%           Functions should be of the form K(x, tau) OR  
%				{K(x, tau, ...), varargin{:}}, etc.
%
% tau -> The vector of your output function (e.g. T1, T2, ec) Usually
%			logspace. Size [m x 1]
%				
% Outputs:
% out -> Struct containing the inverted spectrum.
%   out.t = First dimension in transformed (K) space. [size mx1]
%   out.f = Spectra in the transformed space [size mx1]
%   out.F = The spectrum in the compressed domain.
%   out.K = The kernel functions (combined) in the compressed domain.
%   out.alpha = The list of alphas used.
%
% Usage:
% out = linvert(data, K[, options_struct, ... (options)]);

d = data;

if(~isempty(varargin) && isstruct(varargin{1}))
	% In this case we'll assume you've given us a valid opts struct.
	o = varargin{1};
	
	if(length(varargin) > 1)
		o = laplaceOptions(o, varargin{2:end});
	end
else
	o = laplaceOptions(varargin{:});
end

% Parse the options out into individual variables:
svdd = length(d.x);
if(o.nSVDs <= 0)
	o.nSVDs = svdd;
end

nSVDs = o.nSVDs;
nAlphas = o.nAlphas;
alpha = o.alpha;
a_log = o.alpha_log;
a_end = o.alpha_end;
a_mode = o.alpha_mode;
a_list = o.alpha_list;
alpha_conv = o.alpha_conv;
verbose = o.verbose;
xscale = o.xscale;
yscale = o.yscale;
xlab = o.xlabel;
ylab = o.ylabel;
dataperc = o.dataperc;
opts = o;
o = o.optims;

% Generate the matrices for the kernel functions.
if(iscell(K))
	k = K{1};
	k = k(d.x, tau, K{2:end});
else
	k = K(d.x, tau);
end

% Get the SVDs now:
[U, S, D] = svds(k, nSVDs(1));

% The "S" matrix is the diagonal matrix of singular values, organized in
% descending order, and are a representation of the magnitude of each
% eigenvalue in the Kernel vector. We can eliminate all values of S which
% are below a certain fraction.
S = diag(S);
if dataperc > 0
	con = S > S(1)*dataperc;
else
% 	Mc = U'*d.y;
% 	c = logical(abs(Mc) >= d.std);
% 	con = c;
	con = S.^2 > (S(1)/d.std).^2;
end

n = sum(con(:));

% Smooth truncation.

S = S(con);
[S, i] = sort(S, 'descend'); % Sort this like a proper SVD

% Get D in the compressed subspace.
D = D(:, con);
D = D(:, i);   % Sort 

% U is the unitary transforms which transform between the SVD and
% data bases. Use an inverse transform on the data to bring the data into
% the SVD space, then select the subset of the data which are in our
% compressed space.
Mc = U'*d.y;
Mc = Mc(con);
Mc = Mc(i);		% MCI was a a telecommunications company before it changed 
					% its name to WorldCom. It is now owned by Verison.

% Create the kernel function in the compressed space by recombining S and
% D. S can be constructed this way because the singular value matrix is
% always a diagonal matrix with the singular values along the diagonal.
K = diag(S)*D';
id = eye(n, n);	% An identity matrix.

% Set up the output structure.
out.t = tau;		
out.f = {};
out.c = {};
out.K = K;
out.kf = k;
out.alpha = [];
out.alpha_opt = [];
out.conv = 0;
out.opts = opts;
out.ds = d;

st = n*(d.std).^2;

% Now set up the guess vector.
a = find(strcmp(varargin, 'guess'), 1, 'first');
C0 = ones(n, 1); % Initial guess

if(~isempty(a))
	f = varargin{a+1};
	sf = size(f);

	
	if(isnumeric(f) && (all(sf == ps) || all(sf == fliplr(ps))))
		if(sf ~= ps)
			f = f';
		end
		
		C0 = -(f - Mc)/alpha;
		o.TypicalX = C0;
	else
		warning('The best guess provided was not valid, using default.'); %#ok
	end
else
	% Calculate the first best guess using the first alpha.
	C0 = fminunc(@(c)minvert(c, Mc, K, alpha, id), C0, o);
end

conv = 0;

figname = 'Laplace Inversion';
h = findobj('type', 'figure', 'name', figname);
if(isempty(h))
	h = figure('name', figname);
end

% Reset window size.
set(h, 'Units', 'normalized');
hpos = get(h, 'Position');
wh = 0.65;					% Window height, normalized units
ww = 0.4;					% Window width, normalized units

% Colors
bc = [hex2dec('39'), hex2dec('53'), hex2dec('a4')];
rc = [hex2dec('b5'), hex2dec('13'), hex2dec('0c')];
gc = [hex2dec('00'), hex2dec('79'), hex2dec('35')];

ca = caxis;
cscale = (max(ca)-min(ca))/hex2dec('ff');
coff = min(ca);

bc = bc*cscale + coff;
rc = rc*cscale + coff;
gc = gc*cscale + coff;

hpos(2:end) = [hpos(2)+hpos(4)-wh, ww, wh];
set(h, 'Position', hpos);

brd = false;			% Right now a single boolean is available, but this is
lcurve = false;		% kept separate for future-proofing when additional 
							% methods of regularization parameter choice are
							% added.

if strcmp(a_mode, 'BRD')
	brd = true;
elseif strcmp(a_mode, 'lcurve')
	if isempty(a_list)
		if a_log
			a_list = logspace(alpha, a_end, nAlphas);
		else
			a_list = linspace(alpha, a_end, nAlphas);
		end
	end
	
	lcurve = true;
	
	nAlphas = length(a_list);
else
	error('Invalid alpha mode.');
end

out.eta = zeros(1, nAlphas);
out.rho = zeros(1, nAlphas);

for i = 1:nAlphas
	if brd
		out.alpha(i) = alpha;
	elseif lcurve
		alpha = a_list(i);
		out.alpha(i) = alpha;		% Might be better to pre-allocate, but 
											% I've put this here in case a future
											% version of this can 'fail-safe' and
											% return an interrupted output.
	else
		error('How did you even get this far with an invalid alpha mode?')
	end
	c = C0;
	eflag = -1;
	
	t = cputime;
% 	while eflag < 1
% 		[c, ~, eflag] = fminunc(@(c)minvert(c, Mc, K, alpha, id), c, o);
% 	end

	c = tikhonov(U, S, D, c, alpha);

	t2 = cputime-t;
	minute = floor(t2/60);
	second = t2-minute*60;
		
	C0 = c;
	f = max(0, K'*c);
	
	out.c{i} = c;
	out.f{i} = f;
	
	if(verbose)
		fprintf('%d - Elapsed time %02.0g minutes %02.0g seconds.\n', ...
			i, minute, second);
	end
	
	h2 = gcf;

	set(0, 'CurrentFigure', h);
	subplot(4, 1, [1 2]);
	
	plot(out.t, out.f{i}, 'LineWidth', 1, 'Color', bc, ...
		'LineSmoothing', 'on');
	
	set(gca, 'XScale', xscale);
	set(gca, 'YScale', yscale);
	
	title(sprintf('Not yet converged at %d: \\alpha = %02.2g', i, alpha));
	xlabel(xlab);
	ylabel(ylab);
	
	subplot(4, 1, 3);
	eta = log(norm(f, 'fro'));
	rho = log(alpha*norm(c, 'fro'));

	out.eta(i) = eta;
	out.rho(i) = rho;

	plot(out.eta(1:i), out.rho(1:i), 'o:', ...
		'LineWidth', 2, 'MarkerEdgeColor', 'k', ...
		'MarkerSize', 10, 'MarkerFaceColor', 'w', 'Color', rc, ...
		'LineSmoothing', 'on');
	xlabel('\eta');
	ylabel('\rho');	
		
	subplot(4, 1, 4);	
	if(i >= 4)
		out.rho2d = derivative(out.eta(1:i), out.rho(1:i), 2, 1, true);
		semilogx(out.alpha(1:i), out.rho2d, 'gx--', 'LineWidth', 2, ...
			'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'Color', gc, ...
			'LineSmoothing', 'on');
		
		xlabel('\alpha')
		ylabel('d\rho^2/d^2\eta');
	end
	
	drawnow;

	set(0, 'CurrentFigure', h2);
	
	if brd
		% Compute the new optimal alpha.
		a_opt = sqrt(st)/norm(c, 'fro');

		out.alpha_opt(i) = a_opt;
		if(abs((alpha-a_opt)/alpha) < alpha_conv) % Convergence condition.
			conv = 1;
			out.conv = 1;

			title(sprintf('Converged at %d: \\alpha = %02.2g', i, alpha));
			break;
		end

		alpha = a_opt;
	end
	
	
end

if i < nAlphas
	out.eta = out.eta(1:i);
	out.rho = out.rho(1:i);
end

if(length(out.rho) >= 4)
	out.rho2d = derivative(out.eta, out.rho, 2, 1, true);
else
	out.rho2d = [];
end

if(verbose)
	if(conv)
		fprintf('alpha converged at %3.3f after %d attempts.\n', alpha, length(out.alpha));
	else
		fprintf('alpha failed to converge after %d attempts.\n', nAlphas);
	end
end