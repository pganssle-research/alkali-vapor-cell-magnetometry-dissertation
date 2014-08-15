%% Generate Data
t = linspace(0, 12, 30);
V = linspace(0, 0.38, 29);

% T2 Decay
% Decane T2 - 1.32 s
% Water T2 - 2.47 s
dt2 = 1.32;
wt2 = 2.47;

data_t2 = (0.4*exp_decay(t, dt2) + 0.6*exp_decay(t, wt2))';

% Diffusion decay
% Decane D = 1.80 (40 C)
% Water D = 3.366 (40 C)
dD = 1.8;
wD = 3.366;
data_v = (0.4*diff_kernel(V, dD, 6, 0.075, 4257) + 0.6*diff_kernel(V,wD, 6, 0.075, 4257))';

noise_st = 3;
data_amp = 250;

data = data_amp*data_t2*data_v;

data_t2n = data_amp*data_t2 + normrnd(0, noise_st, size(data_t2));
data_Dn = data_amp*data_v + normrnd(0, noise_st, size(data_v));

dst = make_data_struct(t, data_t2n, noise_st);
dsd = make_data_struct(V, data_Dn, noise_st);

%% Set up inversion
opts1dT2 = laplaceOptions('dataperc', 1e-9, 'verbose', false, ...
	'alpha_mode', 'lcurve', 'alpha', 5, 'alpha_end', -8, ...
'xlabel', 'T_2 (s)', 'ylabel', 'Intensity (pT)', ...
'alpha_log', true);

opts1dD = laplaceOptions('dataperc', 1e-9, 'verbose', false, ...
	'alpha_mode', 'lcurve', 'alpha', 5, 'alpha_end', -8, ...
'xlabel', 'Diffusion Coefficient (10^{-5} cm^2 s^{-1}', ...
'ylabel', 'Intensity (pT)', 'alpha_log', true);

K1 = @exp_decay;
K2 = {@diff_kernel, 6, 0.075, 4257}; % {@func, args} - n, tau, gamma

tau1lin = linspace(0.5, 5, 100);
tau2lin = linspace(0.5, 5, 100);

tau1log = logspace(log10(0.5), log10(5), 100);
tau2log = logspace(log10(0.5), log10(5), 100);


%% Run whichever inversions you want:

run_1d_t2_lin = true;
run_1d_t2_log = false;
run_1d_D_lin = false;
run_1d_D_log = false;

%%
if run_1d_t2_lin
	%% Putting this in a cell so you can run it independently
	% 1D inversion, t2, linear tau. - Remove the );% from this line and
	% uncomment the next line if you want to change the alpha sweep.
	o1tli = linvert(dst, K1, tau1lin, opts1dT2);%, ...
% 		'alpha', -4, 'alpha_end', -10);
	
	%%
end

if run_1d_t2_log
	%% Putting this in a cell so you can run it independently
	% 1D inversion, t2, linear tau. - Remove the );% from this line and
	% uncomment the next line if you want to change the alpha sweep.
	o1tli = linvert(dst, K1, tau1log, opts1dT2);%, ...
		%'alpha', 2, 'alpha_end', -7);
	
	%%
end

if run_1d_D_lin
	%% Putting this in a cell so you can run it independently
	% 1D inversion, t2, linear tau. - Remove the );% from this line and
	% uncomment the next line if you want to change the alpha sweep.
	o1Dli = linvert(dsd, K2, tau2lin, opts1dD);%, ...
		%'alpha', 2, 'alpha_end', -7);
	
	%%
end


if run_1d_t2_lin
	%% Putting this in a cell so you can run it independently
	% 1D inversion, t2, linear tau. - Remove the );% from this line and
	% uncomment the next line if you want to change the alpha sweep.
	o1Dlo = linvert(dsd, K2, tau2log, opts1dD);%, ...
		%'alpha', 2, 'alpha_end', -7);
	
	%%
end