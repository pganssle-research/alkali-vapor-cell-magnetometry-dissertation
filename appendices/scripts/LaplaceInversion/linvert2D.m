function out = linvert2D(data, K1, K2, tau1, tau2, varargin)
% This function takes a data structure (x, y, z) and calculates its 2D
% Laplace inversion with Tikhonov inversion. This is an ill-posed problem,
% and boils down to executing the minimization:
%
%	argmin[F] (||Z - K1FK2'||^2 + a*||F||)
%
% As this could be a VERY time-consuming process with large data-sets, the
% data are first compressed to include only the singular values of the
% kernel functions with the largest probability density.
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
%   data.y = Second dimension (e.g. Gradient strength or time) [size mx1]
%   data.z = Measurement response [size nxm]
%
% K1, K2 -> Kernel functions. Can be anonymous function or function handles
%           Functions should be of the form K1(x, tau1), K2(y, tau2) OR
%           {K1(x, tau1, ...), varargin{:}}, etc.
%
% tau1, tau2 -> The vectors in the transformed space. As the transform is
%						performed in a compressed data space anyway, these can be
%						very large without causing problems. According to the
%						paper, they 
%
% Outputs:
% out -> Struct containing the inverted spectrum
%   out.t1 = First dimension in transformed (K1) space. [size Nx1]
%   out.t2 = Second dimension in transformed (K2) space. [size Mx1]
%   out.f = Spectrum in the transrmed space [size NxM]
%   out.F = The spectrum in the compressed domain.
%   out.K = The kernel functions (combined) in the compressed domain.
%				This is useful for if you want to recreate the spectrum later
%				(for guesses and such).
%
% Usage:
% out = processData(data, K1, K2, ... (options));

d = data;

if ~isempty(varargin) && isstruct(varargin{1})
	% In this case we'll assume you've given us a valid opts struct.
	o = varargin{1};
	
	if(length(varargin) > 1)
		o = laplaceOptions(o, varargin{2:end});
	end
else
	o = laplaceOptions(varargin{:});
end

% Parse the options out into individual variables:
svdd = [length(d.x), length(d.y)];
if(length(o.nSVDs) < 2)
	o.nSVDs(2) = svdd(2);
end
o.nSVDs(o.nSVDs <= 0) = svdd(o.nSVDs <= 0);

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
if(iscell(K1))
	k1 = K1{1};
	k1 = k1(d.x, tau1, K1{2:end});
else
	k1 = K1(d.x, tau1);
end

if(iscell(K2))
	k2 = K2{1};
	k2 = k2(d.y, tau2, K2{2:end});
else
	k2 = K2(d.y, tau2);
end

% Get the SVDs now:
[U1, S1, D1] = svds(k1, nSVDs(1));
[U2, S2, D2] = svds(k2, nSVDs(2));

% The "S" matrix is the diagonal matrix of singular values, organized in
% descending order, and are a representation of the magnitude of each
% eigenvalue in the Kernel vector. We can eliminate all combinations of S1
% and S2 which are below a certain fraction of the largest singular value
% (S1S2). These will be the singular values in our compressed space.
sd = diag(S1)*diag(S2)';

if dataperc > 0
	con = sd > S1(1, 1)*S2(1, 1)*dataperc;
else
	Mc = U1'*d.z*U2;
	c = logical(abs(Mc) >= d.std);
	con = c;
end

n = sum(con(:));

S = sd(con);
[S, i] = sort(S, 'descend'); % Sort this like a proper SVD

% We want to vectorize our data matrix so that the calculations can go
% smoothly. We'll use the equality vec(K1*F*K2') = kron(K2, K1)*vec(F) to
% preserve the operations of the kernel in the vectorized space.
D = kron(D2, D1);
D = D(:, con);
D = D(:, i);   % Sort 

% U1 and U2 are the unitary transforms which transform between the SVD and
% data bases. Use an inverse transform on the data to bring the data into
% the SVD space, then select the subset of the data which are in our
% compressed space.
Mc = U1'*d.z*U2;
Mc = Mc(con);
Mc = Mc(i);		% MCI was a a telecommunications company before it changed 
					% its name to WorldCom. It is now owned by Verison.

% Create the kernel function in the compressed space by recombining S and
% D. S can be constructed this way because the singular value matrix is
% always a diagonal matrix with the singular values along the diagonal.
K = diag(S)*D';
id = eye(n, n);	% An identity matrix.

svd1 = struct('U', U1, 'S', S1, 'D', D1);
svd2 = struct('U', U2, 'S', S2, 'D', D2);

% Set up the output structure.
out.t1 = tau1;		
out.t2 = tau2;
out.f = {};
out.c = {};
out.K = K;
out.kf1 = k1;
out.kf2 = k2;
out.conv = 0;
out.opts = opts;
out.ds = d;
out.svd1 = svd1;
out.svd2 = svd2;
out.alpha = [];

ps = zeros(1, 2);
ps(1) = length(tau1);
ps(2) = length(tau2);

st = n*d.std;

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

		C0 = -(K*reshape(f, prod(ps), 1) - Mc)/alpha;
		o.TypicalX = C0;
	else
		warning('The best guess provided was not valid, using default.'); %#ok
	end
end

ps = num2cell(ps);

% Subplot size
conv = 0;
figname = 'Laplace Inversion';

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
	
	t = cputime;
	eflag = -1;
	c = C0;
	
	while eflag < 1
		[c, ~, eflag] = fminunc(@(c)minvert(c, Mc, K, alpha, id), c, o);
	end
	
	t2 = cputime-t;
	minute = floor(t2/60);
	second = t2-minute*60;

	C0 = c;
	f = max(0, reshape(K'*c, ps{:}));

	out.c{i} = c;
	out.f{i} = f;

	if(verbose)
		fprintf('%d - Elapsed time %02.0g minutes %02.0g seconds.\n', i, minute, second);
	end
	
	h2 = gcf;
	h = findobj('type', 'figure', 'name', figname);
	if(isempty(h))
		figure('name', figname);
	else
		set(0, 'CurrentFigure', h);
	end
	
	contour(out.t1, out.t2, out.f{i}');

	set(gca, 'YScale', yscale);
	set(gca, 'XScale', xscale);

	title(sprintf('%d: \\alpha = %02.2g', i, alpha));
	
	xlabel(xlab);
	ylabel(ylab);
	if ishandle(h2)
		set(0, 'CurrentFigure', h2);
	end
	drawnow;

	if brd
		% Compute the new optimal alpha.
		a_opt = sqrt(st)/norm(c, 'fro');
		if(abs((alpha-a_opt)/alpha) < alpha_conv) % Convergence condition.
			conv = 1;
			out.conv = 1;
			break;
		end

		alpha = a_opt;
	end
end

if(verbose)
	if(conv)
		fprintf('alpha converged at %3.3f after %d attempts.', alpha, length(out.alpha));
	else
		fprintf('alpha failed to converge after %d attempts.', nAlphas);
	end
end

if(conv)
	clf;

	contour(out.t1, out.t2, out.f{i}');

	set(gca, 'YScale', yscale);
	set(gca, 'XScale', xscale);

	title(sprintf('\\alpha = %02.2g', alpha));

	xlabel(xlab);
	ylabel(ylab);
end