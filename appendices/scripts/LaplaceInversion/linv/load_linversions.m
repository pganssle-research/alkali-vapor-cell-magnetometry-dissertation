% This worksheet can be changed according to your liking and is largely an
% example of how to get this is done. Make sure that whatever you called 
% '2012-09-21-DataForBoqin.mat' is in your path when you execute it.

%% Load the data
data_file = '2012-09-21-DataForBoqin.mat';
load(data_file, 'x', 'y', 'z', 'stdev');
ds = make_data_struct(x, y, z, stdev);

%% Set the data percentage parameter.
dperc = 1e-3;


%% Prepare the kernel functions
K1 = @exp_decay;
K2 = {@diff_kernel, 6, 0.075, 4257}; % {@func, args} - n, tau, gamma

%% Prepare the LaplaceOptions struct
opts = laplaceOptions('dataperc', dperc, 'verbose', false, ...
	'nAlphas', 100, 'alpha', 1e-6, 'alpha_conv', 0.001,...
	'xlabel', 'T_2 (s)', ...
	'ylabel', 'Diffusion Coefficient (10^{-5} cm^2 s^{-1})');

%% Generate the logarithmic tau spacing
tau1log = logspace(log10(0.5), log10(4), 150);
tau2log = logspace(log10(0.5), log10(4), 150);

%% Run the logarithmic tau spacing
opts.dataperc = dperc;
opts.xscale = 'log';
opts.yscale = 'log';
o2log = linvert2D(ds, K1, K2, tau1log, tau2log, opts);

%% Generate the linear tau spacing
tau1lin = linspace(0.5, 4, 150);
tau2lin = linspace(0.5, 4, 150);

%% Run the linear tau spacing
opts.dataperc = dperc;
opts.xscale = 'linear';
opts.yscale = 'linear';
o2lin = linvert2D(ds, K1, K2, tau1lin, tau2lin, opts);

%% This next section doesn't seem to be working very well - running the 
%% script with truncation at the noise level.

%% Linear with noise truncation.
opts.dataperc = -1;
opts.xscale = 'linear';
opts.yscale = 'linear';
o2lin2 = linvert2D(ds, K1, K2, tau1lin, tau2lin, opts);

%% Logarithmic with noise truncation
opts.dataperc = -1;
opts.xscale = 'log';
opts.yscale = 'log';
o2log2 = linvert2D(ds, K1, K2, tau1log, tau2log, opts);

%% Open the plotter
plot_linversion(o2log);
% plot_linversion(o2log2);
% plot_linversion(o2lin);
% plot_linversion(o2lin2);



