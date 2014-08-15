function opt = laplaceOptions(varargin)
% Function for generating an options structure for use with the
% linvert2D and linvert functions. Generates a structure which can be
% re-used and modified directly.
%
% Options (use 'name', value)
% 'nSVDs' -> Maximum number of SVDs to use in this calculation. Pass either
%				a scalar, which will be used for both kernels, or a 2x1 vector,
%           where nSVDs(1) is number of SVDs for K1, nSVDs(2) is the number
%           of SVDs for K2. Default is [length(data.x), length(data.y)].
%				Pass 0 or -1 per dimension for default.
%
% 'nAlphas' -> Number of alphas to try out Default 36
%
% 'alpha'	-> Initial alpha choice (default 1)
%
% 'alpha_conv' -> Condition for alpha convergence (as a percentage from
%						optimum). Default 0.05
%
% 'xscale' -> Pass either 'log' or 'linear', default: linear
%
% 'yscale' -> Pass either 'log' or 'linear' default: linear
%
% 'xlabel' -> Label for the x axis (no default)
%
% 'ylabel' -> Label for the y axis (no default)
%
% 'verbose' -> Boolean - if you want printouts of progress. Default: true
%
% 'dataperc' -> Lower limit on the ratio of an SVD included in the
%						compression with the highest SVD. Default 1e-4; 
%
% 'guess'	-> An initial guess for the spectrum - should be size
%					[length(tau1), length(tau2)].
%
% 'optims'   -> Pass an optimset that will be (indirectly) passed to the
%					 fmincon. Default values set by this function are:
%
%					'GradObj' = 'on'
%					'Hessian' = 'on'
%					'TolX' = 1e-9
%					'TolFun' = 1e-11
%					'Display', 'off'
%

% Establish the default options.
svdd = [0, 0];
svals = {'linear', 'log'};
ops = optimset('GradObj', 'on', ...
					'Hessian', 'on', ...
					'TolX', 1e-11, ...
					'TolFun',  1e-12, ...
					'LargeScale', 'on',  ...
					'MaxIter', 1000, ...
					'MaxFunEvals', 1000, ...
					'Display', 'off');


opt = struct('nSVDs', svdd, 'nAlphas', 36, 'alpha', 1, ...
	'alpha_conv', 0.05, 'verbose', 1, 'xscale', 'linear', ...
	'yscale', 'linear', 'xlabel', '', 'ylabel', '', ...
	'dataperc', 1e-4, 'optims', ops, 'alpha_range', [], ...
	'alpha_scale', 'lin');


if(isstruct(varargin{1}))
	o = varargin{1};
	varargin(1) = [];
	
	opt = merge_struct(opt, o);
end


% Process options first.
% Max number of SVDs - default is max num possible.
ns = find(strcmp(varargin, 'nSVDs'));
if(~isempty(ns))
	NS = varargin{ns(1)+1};
	
	if(~isnumeric(NS))
		error('Number of SVDs must be numeric.');
	else
		if(isscalar(NS))
			opt.nSVDs(:) = NS;
		else
			opt.nSVDs(:) = NS(1:2);
		end
	end

	varargin([ns, ns+1]) = [];
end

% Number of alphas
a = find(strcmp(varargin, 'nAlphas'));
if(~isempty(a))
	na = varargin{a(1)+1};
	if(isnumeric(na) && isscalar(na) && na > 0)
		opt.nAlphas = na;
	else
		warning('Invalid nAlphas - using default value of %d.', ...
					opt.nAlphas); %#ok
	end
	
	varargin([a, a+1]) = [];
end

% Initial Alpha
a = find(strcmp(varargin, 'alpha'), 1, 'first');
if(~isempty(a))
	al = varargin{a(1)+1};
	if(isscalar(al) && isnumeric(al) && al >= 0)
		opt.alpha = al;
	else
		warning('Invalid alpha guess, using default value of %d.', ...
					opt.alpha); %#ok
	end
	
	varargin([a, a+1]) = [];
end

% Alpha convergence level
a = find(strcmp(varargin, 'alpha_conv'));
if(~isempty(a))
	na = varargin{a(1)+1};
	if(isnumeric(na) && isscalar(na) && na > 0)
		opt.alpha_conv = na;
	else
		warning('Invalid alpha convergence, using default value of %d.',...
						opt.alpha_conv); %#ok
	end
	
	varargin([a, a+1]) = [];
end

% Verbosity
a = find(strcmp(varargin, 'verbose'));
if(~isempty(a))
	v = varargin{a(1)+1};
	if(~isempty(v))
		opt.verbose = logical(v);
	else
		if(opt.verbose)
			str = 'true';
		else
			str = 'false';
		end
		
		warning('Invalid verbosity, using default value of ''%s''.',  ...
					str);		%#ok
	end
	
	varargin([a, a+1]) = [];
end

% XScale
a = find(strcmp(varargin, 'xscale'));

if(~isempty(a))
	v = varargin{a(1)+1};
	if(~isempty(v))
		v = svals{strcmp(v, svals)};
	end
	
	if(~isempty(v))
		opt.xscale = v;
	else
		warning('Invalid xscale, using default value of ''%s''', ...
					opt.xscale); %#ok
	end
	
	varargin([a, a+1]) = [];		
end

% YScale
a = find(strcmp(varargin, 'yscale'));

if(~isempty(a))
	v = varargin{a(1)+1};
	if(~isempty(v))
		v = svals{strcmp(v, svals)};
	end
	
	if(~isempty(v))
		opt.yscale = v;
	else
		warning('Invalid yscale, using default value of ''%s''', ...
			opt.yscale); %#ok
	end
	
	varargin([a, a+1]) = [];		
end

% XLabel
a = find(strcmp(varargin, 'xlabel'));

if(~isempty(a))
	v = varargin{a(1)+1};
	
	if(ischar(v))
		opt.xlabel = v;
	end
	
	varargin([a, a+1]) = [];
end

% YLabel
a = find(strcmp(varargin, 'ylabel'));

if(~isempty(a))
	v = varargin{a(1)+1};
	
	if(ischar(v))
		opt.ylabel = v;
	end
	
	varargin([a, a+1]) = [];
end

% Data percentage
a = find(strcmp(varargin, 'dataperc'));
if(~isempty(a))
	dp = varargin{a(1)+1};
	if(~isempty(dp) && isnumeric(dp))
		opt.dataperc = dp(1);
	else
		warning('Invalid data percentage, using default value of %d', ...
			opt.dataperc); %#ok
	end
	
	varargin([a, a+1]) = [];	
end

% Optimization parameters
a = find(strcmp(varargin, 'optims'));
if(~isempty(a))
	os = varargin{a(1)+1};
	
	if(~isempty(os) && isstruct(os))
		o = os;
	else
		warning('Invalid optimset, using default values.'); %#ok
	end
end

a = find(strcmp(varargin, 'alpha_range'));
if(~isempty(a))
	os = varargin{a(1)+1};
	if(~isempty(os) && isvector(os))
		opt.alpha_range = os;
	else
		warning('Invalid alpha range, using default values.'); %#ok
	end
end

a = find(strcmp(varargin, 'alpha_scale'));
if(~isempty(a))
	os = varargin{a(1)+1};
	if(~isempty(os) && ischar(os))
		opt.alpha_scale = os;
	else
		warning('Invalid alpha scale, using default value of ''lin''.'); %#ok
	end
end


% Set default values where unset.
if(exist('o', 'var'))
	o = merge_struct(ops, o);
else
	o = ops;
end

opt.optims = o;
