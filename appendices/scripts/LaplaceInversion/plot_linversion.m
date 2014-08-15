function us = plot_linversion(linv_spec)
% Function for plotting the results of laplace inversions, whether or not
% they have converged, and navigating the data.
%
% Usage:
% plot_linversion(linv_spec);

% Make or locate the figure
figname = 'Laplace Inversion Plotter';
fig = findobj('type', 'figure', 'name', figname);
if(isempty(fig))
	fig = figure('name', figname);
elseif(length(fig) > 1)
	fig = fig(1);
end

if iscell(linv_spec)
	s = linv_spec{1};
	cell_spec = linv_spec;
else
	s = linv_spec;
	cell_spec = [];
end

s.n = length(s.alpha);
s = calc_fits_and_cost(s);
s = calc_lcurve(s);

us = build_ui(s, fig, cell_spec);

set(us.vpull, 'value', s.n);
change_num(0, 0, fig);

function the_title = get_title(s)
if s.conv
	the_title = sprintf(...
	'Laplace Inversion Converged (%d steps). Step %%d: \\\\alpha = %%02.2g', ...
	s.n);
else
	the_title = sprintf(...
	'Inversion failed to converge (%d steps). Step %%d: \\\\alpha = %%02.2g', ...
	s.n);
end

function s = calc_fits_and_cost(s)
% Calculate how the fits look.
if(isfield(s.ds, 'z'))
	s.fits = cellfun(@(x)s.kf1*x*s.kf2', s.f, 'UniformOutput', false);
	s.cost = calc_cost(s.ds.z, kron(s.kf2, s.kf1), s.f);
else
	s.fits = cellfun(@(x)x'*s.kf', s.f, 'UniformOutput', false);
	s.cost = calc_cost(s.ds.y, s.kf, s.f);
end

function change_data(~, ~, fig)
% Change what we're looking at
%
% 1 = Spectrum
% 2 = Data
% 3 = Alphas
% 4 = Alpha Difference
% 5 = L-Curve
us = get(fig, 'UserData');
val = get(us.data_mode, 'value');
strs = get(us.disp_mode, 'string');

nstrs = '';
for i = 1:size(strs, 1)
	nstrs = [nstrs, deblank(strs(i, :)), '|']; %#ok;
end
nstrs(end) = [];
strs = nstrs;

set(fig, 'UserData', us);
us = arrange_pulldowns(fig);

if val == 1
	if(~strcmp(strs, us.dmstrs{1}))
		type = get(us.disp_mode, 'value');
		if(type > sum(us.dmstrs{1} == 1)+1)
			set(us.disp_mode, 'value', 1);
		end
		
		set(us.disp_mode, 'string', us.dmstrs{1});
	end
	change_num(0, 0, fig);
elseif val == 2
	if(~strcmp(strs, us.dmstrs{2}))
		type = get(us.disp_mode, 'value');
		if(type > sum(us.dmstrs{2} == 1)+1)
			set(us.disp_mode, 'value', 1);
		end
		
		set(us.disp_mode, 'string', us.dmstrs{2});
	end
	
	change_num(0, 0, fig);
elseif val >= 3 && val <= 5
	% This is where we plot the actual alphas.
	% Check the strings
	if(~strcmp(strs, us.dmstrs{3}))
		set(us.disp_mode, 'string', us.dmstrs{3});
		set(us.disp_mode, 'value', 1);
		
		if(val == 3)
			set(us.vpull, 'value', 2);
		else
			set(us.vpull, 'value', 1);
		end
	end
	
	st = get(us.vpull, 'value');
	en = us.s.n;
	
	if val == 3 || val == 4
		type = us.alphas_type;
	else
		type = us.cost_type;
	end
	
	if val == 3
		data = us.s.alpha(st:en);
		x = 1:length(data);
	elseif val == 4
		data = 100*(us.s.alpha-us.s.alpha_opt)./us.s.alpha;
		data = data(st:end);
		x = 1:length(data);
	else
		x = us.s.eta(st:en);
		data = us.s.rho(st:en);
		
		%x = us.s.alpha(st:en);
		%data = us.s.cost(st:en);
	end
	
	set(us.disp_mode, 'Value', type);
	vargs = {'-ok', 'markers', 8, 'MarkerFaceColor', 'k'};
	if(type == 1)		% Lin-Lin
		plot(us.axes, x, data, vargs{:});
	elseif type == 2	% Lin-Log
		semilogy(us.axes, x,  data, vargs{:});
	elseif type == 3 % Log-Lin
		semilogx(us.axes, x,  data, vargs{:});
	elseif type == 4 % Log-Log
		loglog(us.axes, x,  data, vargs{:});
	end
	
	if(val == 3)
		us.title = title(sprintf(...
			'Alphas %d to %d', st, en), ...
			'FontWeight', 'demi');
	elseif val == 4
		us.title = title(sprintf(...
			'Alphas (percent difference from optimal), %d to %d', ...
			st, length(data)),...
		'FontWeight', 'demi');
	elseif val == 5
		us.title = title(sprintf(...
			'L-Curve, %d to %d', ...
			st, en));
	end
end

set(fig, 'UserData', us);

function us = change_type(~, ~, fig)
% Callback for when to change the type (contour/image, etc)
us = get(fig, 'UserData');
val = get(us.data_mode, 'Value');
type = get(us.disp_mode, 'Value');
if val > 2
	if val == 3 || val == 4
		us.alphas_type = type;
	else
		us.cost_type = type;
	end
end

set(fig, 'UserData', us);

if val > 2
	change_data(0, 0, fig);
else
	change_num(0, 0, fig);
end

function change_num(~, ~, fig)
% Callback for when to change the data set

% Determine what we need to plot.
% 1 = Spectrum
% 2 = Data
% 3 = Alphas
us = get(fig, 'UserData');
num = get(us.vpull, 'Value');
val = get(us.data_mode, 'Value');

s = us.s;

if(val == 1 || val == 2)
	type = get(us.disp_mode, 'Value');
	if(val == 1)
		if(us.is2D)
			% Determine the type
			% 1 = Contour
			% 2 = Image

			if(type == 1)
				contour_plot(us, num);
			else
				image_plot(us, num);
			end

			set_labels(s);
		else
			plot(us.axes, s.t, s.f{num});
		end
	elseif(val == 2)
		if(us.is2D)
			% Data
			%
			% Determine the type
			% 1 = Contour
			% 2 = Image
			% 3 = Plots - Dimension 1
			% 4 = Plots - Dimension 2
			d = s.ds;
			if type == 1
				contour(us.axes, d.x, d.y, d.z');
			elseif type == 2
				image(d.x, d.y, d.z, 'CDataMapping', 'scaled');
			elseif type == 3 || type == 4
				data = d.z;
				fit = s.fits{num};
				if type == 3
					x = d.x;
				else
					x = d.y;
				end
				
				plot(us.axes, x, data, '.');
				hold on
				plot(us.axes, x, fit, '-');
				hold off
			end
		else
			plot(us.axes, s.ds.x, s.ds.y, 'ob');
			
			% Also plot the fit.
			hold on
			plot(us.axes, s.ds.x, s.fits{num}, '-b');
			hold off
		end			
	end
	
	us.title = title(sprintf(us.the_title, num, s.alpha(num)), 'FontWeight', 'demi');
else
	change_data(0, 0, fig);
end
set(fig, 'UserData', us);

function axes = contour_plot(us, num)
% Do the contour plot
s = us.s;
axes = contour(us.axes, s.t1, s.t2, s.f{num}');

function axes = image_plot(us, num)
% Do the image plotting
s = us.s;
axes = image(s.t1, s.t2, s.f{num}', 'CDataMapping', 'scaled');
set(gca, 'YDir', 'normal');

function set_labels(s)
if(isfield(s.opts, 'xlabel'))
	xlabel(s.opts.xlabel);
end

if(isfield(s.opts, 'ylabel'))
	ylabel(s.opts.ylabel);
end

function change_cell(~, ~, fig)
% Function for changing which cell you're in
us = get(fig, 'UserData');
val = get(us.cpull, 'Value'); % Which index we're on.
us.cell_spec{us.ind} = us.s;

s = us.cell_spec{val};
prev_ind = us.ind;
us.ind = val;

if ~isfield(s, 'cost')
	s = calc_fits_and_cost(s);
end

if ~isfield(s,'lcurve')
	s = calc_lcurve(s);
end

if ~isfield(s, 'n')
	s.n = length(s.alpha);
end

us.the_title = get_title(s);
us.s = s;

spec_fig = get(us.vpull, 'Value');

set(fig, 'UserData', us);
us = build_ind(fig);

dmode = get(us.data_mode, 'Value');

if us.min_cost && dmode < 3
	[~, i] = min(s.cost);
	set(us.vpull, 'Value', i); 
else
	if spec_fig > s.n || spec_fig == us.cell_spec{prev_ind}.n
		set(us.vpull, 'Value', s.n);
	else
		set(us.vpull, 'Value', spec_fig);
	end
end
update_dimension(fig);
change_num(0, 0, fig);

function ui_struct = build_ui(s, fig, cell_spec)
% Builds the axes and the UI struct telling you which ones are which.
set(0, 'CurrentFigure', fig);
clf;

us = struct('fig', fig, ...	% The figure itself
	'axes', [], ...				% The plot axes
	'title', [], ...				% Plot title
	'vpull', [], ...				% Values pulldown
	'close', [], ...				% Close control
	'legend', [], ...				% The legend
	'legon', [], ...				% Whether or not the legend is on
	'nup', [], ...					% The number to plot at once
	'cax', [], ...					% Current axis
	'data_mode', [], ...			% What data to display.
	'disp_mode', [], ...			% Display mode
	's', [], ...					% The structure
	'dmstrs', [], ...				% Strings for data modes
	'is2D', [], ...				% If it's 2D or 1D.
	'the_title', [], ...			% The base title.
	'cell_spec', [], ...			% The spectrum cells
	'cpull', [], ...			% Pulldown for cell navigation
	'ind', 1, ...
	'scrsz', [], ...
	'pos_vec', [], ...
	'upos_vec', [], ...
	'ups_vec', [], ...
	'cost_type', 4, ...
	'alphas_type', 4, ...
	'min_cost', true ...
	);

the_title = get_title(s);
us.dmstrs = {'Contour|Image', 'Contour|Image|Plots (1)|Plots (2)', ...
	'Lin-Lin|Log-Lin|Lin-Log|Log-Log'};
us.s = s;
us.the_title = the_title;
scrsz = get(0, 'ScreenSize');
figsz = [scrsz(3)*0.125, scrsz(4)*0.125, scrsz(3)*0.6, scrsz(4)*0.6];
set(fig, 'Position', figsz, 'NumberTitle', 'off');

us.axes = axes;
set(us.axes, 'Position', [7/128, 5/64, 1-5/32, 1-1/8]);
us.title = title(sprintf(the_title, 0, 0.0));
set(us.axes, 'Title', us.title);

if ~exist('cell_spec', 'var')
	cell_spec = [];
end

pull_width = 3/32;
pull_height = 0.05;
pull_left = (1-5/32+23/512+1/64);
pull_bot = 0.895;
pull_sep = 0.05;

us.pos_vec = [pull_left, pull_bot, pull_width, pull_height];
us.pull_sep = pull_sep;
us.figsz = figsz;
us.scrsz = scrsz;

% Left, Bottom, Width, Height
pos_vec = us.pos_vec .* [figsz(3), figsz(4), figsz(3), figsz(4)];
ps_vec = [0, pull_sep*figsz(4), 0, 0];

us.upos_vec = pos_vec;
us.ups_vec = ps_vec;

damodepos = pos_vec;
dmodepos = pos_vec - ps_vec;

if ~isempty(cell_spec)
	lspec = length([cell_spec{:}]);
	
	if(isvector(cell_spec))
		format_str = sprintf('%%0%dd|', floor(log10(lspec)));
		inds = num2cel(1:lspec);
	else
		ndig = floor(log10(max(size(cell_spec)))+1);
		format_str = repmat([sprintf('%%0%dd', ndig), ','],...
									1, ndims(cell_spec));
		format_str(end) = '|';
	
		inds = cell(size(size(cell_spec)));
		[inds{:}] = ind2sub(size(cell_spec), 1:lspec);
		inds = num2cell(cell2mat(inds')); % This is a crazy way to do this.		
	end
	
	format_str = repmat(format_str, 1, lspec);
	format_str(end) = [];
	
	cstr = sprintf(format_str, inds{:});
	cpos = pos_vec;
	cpos(2) = (2/32)*figsz(4);
	
	us.cpull = uicontrol('Style', 'popup', 'String', cstr, ...
		'Position', cpos, 'Callback', {@change_cell, fig});
	
	us.cell_spec = cell_spec;
end

dastr = 'Spectrum|Data|Alphas|AlphaDiff|LCurve';
us.data_mode = uicontrol('Style', 'popup', 'String', dastr, ...
	'Position', damodepos);

typestr = us.dmstrs{1};
us.disp_mode = uicontrol('Style', 'popup', 'String', typestr, ...
	'Position', dmodepos);

set(us.disp_mode, 'Callback', {@change_type, fig});
set(us.data_mode, 'Callback', {@change_data, fig});
set(fig, 'UserData', us);
build_ind(fig);
update_dimension(fig);

datacursormode on;

ui_struct = us;

function us = build_ind(fig)
% Call this after the rest of the figure has been constructed to add the
% data structure and whatnot to the UI.
us = get(fig, 'UserData');

pos_vec = us.upos_vec;
ps_vec = us.ups_vec;

if ~isempty(us.vpull)
	delete(us.vpull);
end

tpullpos = pos_vec - 2*ps_vec;
s = us.s;

format_str = sprintf('%%0%dd|', floor(log10(s.n)+1));
format_str = repmat(format_str, 1, s.n);
format_str(end) = [];
inds = num2cell(1:s.n);
tstr = sprintf(format_str, inds{:});

us.vpull = uicontrol('Style', ...
	'popup', 'String', tstr, ...
	'Position', tpullpos);

set(us.vpull, 'Callback', {@change_num, fig});

set(fig, 'UserData', us);

function us = update_dimension(fig)
us = get(fig, 'UserData');
s = us.s;
if isfield(s.ds, 'z')
	us.is2D = true;
else
	us.is2D = false;
end

set(fig, 'UserData', us);

us = arrange_pulldowns(fig);

function us = arrange_pulldowns(fig)
us = get(fig, 'UserData');

pos_vec = us.upos_vec;
ps_vec = us.ups_vec;

val = get(us.data_mode, 'Value');

needs_modes = ((val >= 3) || us.is2D);

if needs_modes
	set(us.disp_mode, 'Position', pos_vec-ps_vec);
	set(us.vpull, 'Position', pos_vec-2*ps_vec);
	set(us.disp_mode, 'Visible', 'on');
else
	set(us.disp_mode, 'Position', pos_vec-2*ps_vec);
	set(us.vpull, 'Position', pos_vec-ps_vec);
	set(us.disp_mode, 'Visible', 'off');
end

set(fig, 'UserData', us);