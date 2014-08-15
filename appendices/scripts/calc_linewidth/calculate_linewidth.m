function [R, RSD, RW, DL] = calculate_linewidth(temp, dims, pps, varargin)
% Linewidth Calculator
% Last Updated: April 2013
% Author: Paul J. Ganssle, Pines Lab, UC Berkeley, Berkeley, CA
%
% This function calculates the fundamental linewidth due to relaxation of
% pumped, uncoated, buffered alkali vapor cells as a function of 
% polarization, temperature and volume. It assumes that the cell is 
% cuboidal. The default values are for a Rb-vapor cell with a combination 
% of neon and nitrogen buffer gas. 
%
% Isotope effects are taken into account only by taking a weighted average
% (weighted by abundance) of the relaxation rates, ignoring any
% second-order mixing effects.
%
% The following works are cited in the comments:
% [1] Killian, T. Physical Review 1926, 27, 578–587
% [2] Wu, Z.; Kitano, M.; Happer, W.; Hou, M.; Daniels, 
%										J. Applied Optics 1986, 25, 4483.
% [3] Alcock, C. B.; Itkin, V. P.; Horrigan, M. K.;
%                     Canadian Metallurgical Quarterly 1984, 23, 309–313.
% [4] Seltzer, S. J., Developments in alkali-metal atomic magnetometry, 
%							   	               Princeton University, 2008.
%
% Inputs:
% ------------------------------------------------------------------------
% temp:		[1xn Double] Cell temperature(s), in degrees C.
% dims:		[1x3 Double] Dimensions of the cell volume, in cm.
% pps:		[mxN Double] partial pressures of each of the N components.
%			Each column in this matrix corresponds to a component, each row
%			corresponds to a different set of partial pressures. Most
%			common application will likely be m = 1, N = 2.
%
% opts:		Pass an lw_opts() struct and/or a combination of label-value
%			pairs setting the options for the calculations. Starting at the
%			4th argument, once a vector 
%			See help	lw_opts() for the default options.
%
%			For N > 2 (where N is the number of gaseous components), must 
%			set the option 'comps' to a cell array of the name of the 
%			components. 
%
%			Use lw_opts('get_components') for a list of components. 
%			To add a components, pass a cell array of add_lw_comp(...)
%			objects to the 'add_comps' option. For adding a single
%			component, a single add_lw_comp() object is acceptable.
% ------------------------------------------------------------------------
%
% Outputs:
% ------------------------------------------------------------------------
% R:		Overall linewidth (in cycles/s a.k.a. Hz)
% RSD:		Linewidth from spin-destruction. (Hz)
% RW:		Linewidth from wall collisions. (Hz)
% DL:		The characteristic diffusion length, in CM.
%
% Note: The cardinality of all the outputs is squeeze(mxnxI) where m is
% the number of rows in pps, n is the number of columns in temp and I is
% the number of isotopes with non-zero abundance in the experiment. The
% ------------------------------------------------------------------------
%
% Usage:
% [R, RSD, RW, DL] = calculate_linewidth(temp, dims, pps[, lwopts]);

if isempty(varargin)
	lwopts = lw_opts();		% Retrieves the default options structure.
else
	lwopts = lw_opts(varargin{:});		% Pass this on.
end

if ~exist('sdcs', 'var') || sdcs == 0
   % By default it's neon and nitrogen
   sdcs = [1.9e-23, 1e-22]; % Spin-destruction cross-sections in cm^2
   
   if length(pps) == 1
       sdcs = sdcs(1);
   elseif length(pps) ~= 2
       error(['Must provide spin-destruction cross sections for more ', ...
			 'than 2 components.']);
   end  
end

if ~exist('ms', 'var') || ms == 0
    ms = [3.35082*1e-26, 4.65173e-26]; % Mass in kg, Neon, N2
    
    if(length(pps) == 1)
        ms = ms(1);
    elseif(length(pps) ~= 2)
      error('Must provide masses for more than 2 components.');
    end
end

Pdflt = 0.5;	% Default value - stored here so it's easier to change.
if ~exist('P', 'var')
	P = Pdflt;
elseif P < 0 || P > 1
	warning(['Invalid polarization value %f provided -', ...
					'using default value of %f'], P, Pdflt);
	P = Pdflt;
end

% Convert temperature to Kelvin scale.
temp = temp + 273.15;

% Rb-Rb cross section
% rn is an approximation of the number density, from:
%  via
% 
rn = (1/temp)*10^(21.866+4.312-4040/temp);  % Rb density in cm^(-3)
rm = 1.41923e-25;                           % Rubidium mass in kg.
kb = 1.38*1e-23;                            % Boltzmann constant, J/K

cell_fill_temp = 293;				% Temperature the cell was filled at in K. 
pps = pps * 133.3;					% Convert pressures from torr to Pa.
n = pps/(kb * cell_fill_temp);   % Number density lost here.
n = n * 1e-6;           % In (number) cm^-3

Mr = (1/rm + 1./ms).^(-1);							% Reduced masses.
vb = sqrt((8*kb*temp)./(pi*Mr)) * 100;			% Avg. atomic velocity in cm/s 
vr = sqrt((16*kb*temp)/(pi*rm)) * 100;

% Standard 'optimal' polarization is 1/2. These are the nuclear slowing-
% down factors. More details in: 
%
q_85 =  (38 + 52 * P^2 + 6 * P^4)/(3 + 10*P^2 + 3 * P^4);
q_87 = (6 + 2*P^2)/(1 + P^2);

% These are relaxation rates - they add linearly.
RRSDB = vr * 9e-18 *rn;
RSDB = sum(vb.*sdcs.*n) + RRSDB;
RSD85 = (1/q_85)*RSDB;
RSD87 = (1/q_87)*RSDB;

RSD = [RSD85, RSD87];		% Spin-destruction Rate is in s^-1. 
									% Separate analysis for 85-Rb and 87-Rb.

D = 0; D1 = 0;					% D =  Diffusion coefficient. (cm^2/s)
									% D1 = Diffusion coefficient for component 1.

% Diffusion scales linearly in p, so D0(P/P0) 
if pps(1) > 0					% Only include diff if it's in the cell. 
    D1 = 0.235 * (temp/305).^(3/2) * (760*133.3/pps(1));    % cm^2/s
    D = D1;
end
			
if pps(2) > 0					% If there's
   D2 = 0.19*(temp/273).^(3/2)*(760*133.3/pps(2));
   if D ~= 0
      D = (1./D1 + 1./D2).^(-1); 
   else
       D = D2;
   end
end

RW = D*(sum((pi./dims).^2));	% Relaxation rate due to wall collisions.
R = RW + RSD;				% Total relaxation rate in Hz.
DL = sqrt(D./R);			% Characteristic diffusion length.
R = R/(2*pi);				% Total relaxation rate in rad Hz.