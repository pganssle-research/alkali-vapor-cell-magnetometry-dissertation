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
alpha_conv = o.alpha_conv;
verbose = o.verbose;
xscale = o.xscale;
yscale = o.yscale;
xlab = o.xlabel;
ylab = o.ylabel;
dataperc = o.dataperc;
alpha_range = o.alpha_range;
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

% Run a test to see where a reasonable cutoff would be:
Mc = U'*d.y;
c = logical(abs(Mc) >= d.std);

% The "S" matrix is the diagonal matrix of singular values, organized in
% descending order, and are a representation of the magnitude of each
% eigenvalue in the Kernel vector. We can eliminate all values of S which
% are below a certain fraction.
S = diag(S);
%con = S > S(1)*dataperc;
con = c;
n = sum(con(:));

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

if(~isempty(alpha_range))
	if(strcmp(opts.alpha_scale, 'log'))
		out.alpha = logspace(alpha_range(1), alpha_range(2), nAlphas);
	else
		out.alpha = linspace(alpha_range(1), alpha_range(2), nAlphas);
	end
	out.alpha_opt = out.alpha;
	out.cost = zeros(size(out.alpha));
	alpha = out.alpha(1);
	conv = 1;
	out.conv = 1;
end

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
end

conv = 0;

figname = 'Laplace Inversion';
h = findobj('type', 'figure', 'name', figname);
if(isempty(h))
	h = figure('name', figname);
end

o.MaxFunEvals = 500;

c_fun = @(c, a)minvert(c, Mc, K, a, id);
f_fun = @(a)max(0, K'*fminunc(@(c)c_fun(c, a), C0, o));

o.Hessian = 'off';
o.GradObj = 'off';
o.Display = 'on';
o.TolFun = 1e-15;
o.LargeScale = 'off';
o.Parallel = 'always';

alpha_opt = fminsearch(@(a)calc_cost(d.y, k, f_fun(10^(a))), alpha, o);
% [X, ~, ~, out] = fmincon(@(alpha)calc_cost(d.y, K, max(0, K'*fminunc(@(c)minvert(c, Mc, K, alpha, id), C0, o)), alpha, [], [], [], [], 1e-20, 0.3);

alpha_opt = 10^(alpha_opt);
out.f = f_fun(alpha_opt);
out.alpha_opt = alpha_opt;

subplot(2, 1, 1);
plot(out.t, out.f);

subplot(2, 1, 2);
plot(out.ds.x, out.ds.y, 'ob', out.ds.x, out.kf*out.f, '-b');


% for i = 1:nAlphas	
% 	t = cputime;
% 	c = fminunc(@(c)minvert(c, Mc, K, alpha, id), C0, o);
% 	
% 	t2 = cputime-t;
% 	minute = floor(t2/60);
% 	second = t2-minute*60;
% 		
% 	C0 = c;
% 	f = max(0, K'*c);
% 	
% 	out.c{i} = c;
% 	out.f{i} = f;
% 	
% 	if(verbose)
% 		fprintf('%d - Elapsed time %02.0g minutes %02.0g seconds.\n', i, minute, second);
% 	end
% 	
% 	h2 = gcf;
% 
% 	set(0, 'CurrentFigure', h);
% 	plot(out.t, out.f{i}, 'LineSmoothing', 'on');
% 	
% 	set(gca, 'XScale', xscale);
% 	set(gca, 'YScale', yscale);
% 	
% 	title(sprintf('Not yet converged at %d: \\alpha = %02.2g', i, alpha));
% 	xlabel(xlab);
% 	ylabel(ylab);
% 	
% 	drawnow;
% 
% 	set(0, 'CurrentFigure', h2);
% 	
% 	if(isempty(alpha_range))
% 		Compute the new optimal alpha.
% 		a_opt = sqrt(st)/norm(c, 'fro');
% 
% 		out.alpha(i) = alpha;
% 		out.alpha_opt(i) = a_opt;
% 		if(abs((alpha-a_opt)/alpha) < alpha_conv) % Convergence condition.
% 			conv = 1;
% 			out.conv = 1;
% 
% 			title(sprintf('Converged at %d: \\alpha = %02.2g', i, alpha));
% 			break;
% 		end	
% 		alpha = a_opt;
% 	elseif(i < nAlphas)
% 		alpha = out.alpha(i+1);
% 	end
% end
% 
% if(verbose)
% 	if(conv)
% 		fprintf('alpha converged at %3.3f after %d attempts.\n', alpha, length(out.alpha));
% 	else
% 		fprintf('alpha failed to converge after %d attempts.\n', nAlphas);
% 	end
% end