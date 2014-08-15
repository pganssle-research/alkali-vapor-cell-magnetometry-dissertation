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
% 'alpha_mode' -> Mode of alpha selection. Options are (currently) 'BRD',
%						'lcurve'. Default is 'lcurve'.
%
% 'alpha' -> Initial alpha choice. Unused if alpha_list is used. Default 1.
%
% 'nAlphas' -> Number of alphas to try out. Default is 20 in lcurve mode,
%					50 in BRD mode.
%
% 'alpha_end' ->	The final value of alpha to try. If no value is passed to 
%						alpha_list and alpha_mode is lcurve, then the alpha list
%						will be generated as:
%							alpha_list = linspace(alpha, alpha_end, nAlphas);
%						for alpha_log = false, otherwise:
%							alpha_list = logspace(alpha, alpha_end, nAlphas);
%						Default is -8 for logscale, 1e-8 for linscale.
%
% 'alpha_list' -> Choose the alphas from a pre-set list. Default choice is
%						to use alpha, alpha_end, nAlphas and alpha_log to
%						generate a list. These are all ignored if a list is
%						passed to this parameter.
%
% 'alpha_conv' -> Condition for alpha convergence (as a percentage from
%						optimum) - only used in 'alpha_mode' = BRD Default 0.05
%
% 'alpha_log'  -> Boolean value. For true, alpha is varied on a log scale.
%						For false, linear. Default is true.
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

opt = struct('nSVDs', svdd, 'nAlphas', 50, 'alpha', 1, ...
	'alpha_conv', 0.05, 'alpha_mode', 'lcurve', 'alpha_end', -8, ...
	'alpha_list', [], 'alpha_log', true, 'verbose', 1, ...
	'xscale', 'linear', 'yscale', 'linear', 'xlabel', '', 'ylabel', '', ...
	'dataperc', 1e-4, 'optims', ops);


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

% Alpha mode.
a = find(strcmp(varargin, 'alpha_mode'));
if(~isempty(a))
	am = varargin{a(1)+1};
	if ischar(am) && find(strcmpi(am, {'BRD', 'lcurve', 'l-curve'}), 1, 'first')
		if strcmpi(am, 'BRD')
			am = 'BRD';
		elseif find(strcmpi(am, {'lcurve', 'l-curve'}), 1, 'first')
			am = 'lcurve';
		end		
		
		opt.alpha_mode = am;		
	else
		warning('Invalid alpha_mode passed as parameter. Using default value instead.'); %#ok
	end
end

a = find(strcmp(varargin, 'alpha_log'));
if(~isempty(a))
	os = varargin{a(1)+1};
	
	try
		opt.alpha_log = logical(os);
	catch %#ok
		warning('Alpha_log parameter could not be converted to logical. Using default value of true instead.') %#ok
		opt.alpha_log = true;
	end
end

brd = false;
lcurve = false;
if(strcmp(opt.alpha_mode, 'BRD'))
	brd = true;
elseif(strcmp(opt.alpha_mode, 'lcurve'))
	lcurve = true;
end

% Number of alphas
a = find(strcmp(varargin, 'nAlphas'));
if(~isempty(a))
	na = varargin{a(1)+1};
	if(isnumeric(na) && isscalar(na) && na > 0)
		opt.nAlphas = na;
	else
		if strcmp(opt.alpha_mode, 'BRD')
			opt.nAlphas = 50;
		else
			opt.nAlphas = 20;
		end
		
		warning('Invalid nAlphas - using default value of %d.', ...
					opt.nAlphas); %#ok
	end
	
	varargin([a, a+1]) = [];
end

% Initial Alpha
a = find(strcmp(varargin, 'alpha'), 1, 'first');
if(~isempty(a))
	al = varargin{a(1)+1};
	if(isscalar(al) && isnumeric(al) && (((brd || (lcurve && ~opt.alpha_log)) && al > 0) || lcurve))
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

a = find(strcmp(varargin, 'alpha_end'));
if(~isempty(a))
	os = varargin{a(1)+1};
	if(~isempty(os) && isscalar(os))
		opt.alpha_end = os;
	else
		if opt.alpha_log
			dflt_val = -8;
		else
			dflt_val = 1e-8;
		end
		
		warning(fprintf('Invalid alpha_end parameter, using default value of %f\n', dflt_val));
		opt.alpha_end = dflt_val;
	end
end

% Alpha list
a = find(strcmp(varargin, 'alpha_list'));
if(~isempty(a))
	os = varargin{a(1)+1};
	if(~isempty(os) && isnumeric(os) && isvector(os))
		opt.alpha_list = os;
	else
		opt.alpha_list = [];
		
		warning('Invalid alpha_list parameter, list will be generated from defined range.'); %#ok
	end
end

% Set default values where unset.
if(exist('o', 'var'))
	o = merge_struct(ops, o);
else
	o = ops;
end

opt.optims = o;
