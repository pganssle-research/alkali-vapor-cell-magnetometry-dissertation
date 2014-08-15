function iso = add_iso(element, number, mass, spin, abundance, ...
									nd_kernel, se)
% Isotope Structure Construction
% Last Updated: April 2013
% Author: Paul J. Ganssle, 
%
% Function for adding new isotopes to the linewidth calculation. This
% basically just generates a structure for each isotope. Repeat once for
% each isotope. Make sure that 'element' is the same for isotopes of the
% same element. Adding a new value with the same element and number
% replaces the previous value.
%
% There is also the optional argument 'se' which specifies spin-exchange
% cross-sections with other alkali metals.
%
% Inputs:
% ------------------------------------------------------------------------
% element:		[String] The name of the element (e.g. Rb, Cs)
% number:		[Integer] Mass number of the isotope (e.g. 85, 87)
% mass:			[Double] Exact mass of the element in au.
% spin:			[Integer] For now, only spins 1/2, 3/2, 5/2 and 7/2 are
%					supported, so pass either 1, 3, 5 or 7. 
% abundance:	[Double] Natural abundance of this isotope.
% nd_kernel:	[Function Handle] A function handle for how the number 
%					density varies as a function of temperature. All functions 
%					should use only element-wise operations.
%
%					Example for Rb:
%					kbat = 1.36259e-22; % Boltzmann's cponstant in cm^3*atm/K
%					lkb = -log10(kbat);
%					nd_kernel = @(T)(1./T)*10.^(-lkb+4.312-4040./T);
%					
%					Alternatively, passing a 2-vector of Virial coefficients 
%					[A, B], the functional defaults to:
%					nd_kernel = @(A, B, T)(1./T).*10.^(lkb + A + B./T);
%					
%					In the Rb example above, [A, B] = [4.312, -4040];
%
% Optional inputs (pass in key-value pairs)
% se:				{{'Element1', 'Element2', ...}, [se1, se2, ... ]}
%					{'Element', se} - Either set a cell array of element names
%					all at once to a vector of spin exchange cross-sections or
%					replace them one at a time. Further entries in the cell
%					array for the same element are ignored.
%					
%					These will merge in with the entire element, and are not
%					unique to the isotope. New information trumps old
%					information, put NaN for placeholder and a negative number
%					revert the current value to NaN.
%
% ------------------------------------------------------------------------
%
% Output:
% ------------------------------------------------------------------------
% iso:			A struct containing all the information.
% ------------------------------------------------------------------------
%
% Prototypes:
% iso = add_iso(element, number, mass, spin, abundance, nd_kernel);

% Input validation
if ~ischar(element)
	error('Element must be a string.');
end

cardin = length(number);

if any(round(number) ~= number) ||...
		any(number <= 0) || ~any(isreal(number))
	error('Mass numbers must be non-zero positive real integers.');
end

% Validate the mass and convert it to kg.
if length(mass) ~= cardin
	error('Must have a mass for each isotope.')
elseif any(~isnumeric(mass)) || any(mass <= 0) || any(~isreal(mass))
	error('Masses must be a positive values in au.');
end

pmkg = 1.660539e-27;		% Close to the proton mass, in kg.
mass = mass*pmkg;

% Validate the spins, then convert them to indexes to the slowing-down
% functionals.
if length(spin) ~= cardin
	error('Number of spins must match number of isotopes.');
elseif any(~isnumeric(spin)) || ...
			any(arrayfun(@(x)~any(x == [1 2 3 4]), spin))
	error(['Each spin must be 1 (S = 1/2), 3 (S = 3/2), ', ...
		'5 (S = 5/2) or 7 (S = 7/2).']);
end

qfs = spin.*2-1;			% Converts 1,3,5,7 to the indexing vector 1,2,3,4.

% Note: Fix this line - I think this spin 1/2 thing is wrong.
QF = {@(P)2, @(P)(6+2*P.^2)./(1+P.^2), ...
	@(P)(38+52*P.^2+6*P.^4)./(3+10*P.^2+3*P.^4), ...
	@(P)(22+70*P.^2+34*P.^4 + 2*P.^6)./(1+7*P.^2+7*P.^4+P.^6)};

if length(abundance) ~= cardin
	error('Must have an abundance for each isotope added.');
elseif ~isnumeric(abundance) || any(abundance < 0) || any(abundance > 1)
	error('Abundance must be a number between 0 and 1.');	
end

if sum(abundance) > 1
	warning('Abundance values sum to %f, which is greater than 1.', ...
		sum(abundance));
end

% Parse what they've entered. There are 4 possibilities - a single function
% handle, a cell array of function handles (length cardinality), an array 
% of virial coefficients or a matrix of virial coefficients of either size 
% (2 x cardinality) or (cardinality x 2). There is an edge case where 
% 2 == cardinality - we'll assume that the virial coefficients are the
% rows.
%
% We're not going to validate whether the functions they pass are actually
% functions of T - they can deal with those errors if they come up.
lkb = -log(1.36259e-22);
virial_func = @(A, T)10.^(lkb + A(1) + A(2)./T);

kernel_error = ['Must provide either 1 or N kernels where N is the ', ...
			'number of isotopes you''re currently adding. You can also ',...
			'add a single set of 3 virial coefficients or an Nx3 ', ...
			'matrix of virial coefficients. See help for more details.'];
if isscalar(nd_kernel)
	% If it's scalar, it can only be a function handle. If it's not, that's
	% an error.
	if isa(nd_kernel, 'function_handle')
		nd_kern = repmat({nd_kernel}, 1, cardin);
	else
		error(kernel_error);
	end	
elseif iscell(nd_kernel)
	% If it's a cell, we'll assume that it's not the virial coefficients,
	% because that's a stupid way to store a matrix or a vector anyway.
	if length(nd_kernel) ~= cardin || ...
		any(~cellfun(@(x)isa(x, 'function_handle'), nd_kernel)) 
		error(kernel_error);
	else
		nd_kern = nd_kernel;
	end
elseif isnumeric(nd_kernel)
	% Now it must be either a matrix or vector of virial coefficients.
	sk = size(nd_kernel);
	if length(sk) > 2 || ~any(size(nd_kernel) == 2)
		error(kernel_error);
	end
	
	s1 = size(nd_kernel, 1);
	s2 = size(nd_kernel, 2);
		
	if s1 ~= s2
		% This is the case where they've inverted it.
		if s1 == 2
			nd_kernel = nd_kernel';			
		end		
		s1 = s2;
	elseif cardin ~= 2
		% If the two are equal but they are not equal to the number of
		% isotopes being added, it's an error.
		error(kernel_error);
	end
	
	if s1 == 1
		nd_kern = repmat({@(T)virial_func(nd_kernel, T)}, 1, cardin);
	elseif s1 == cardin
		nd_kern = arrayfun(@(ii)...
			@(T)virial_func(nd_kernel(ii, :), T), 1:s1, ...
			'UniformOutput', false); 
	else
		error(kernel_error);
	end		
else
	error(kernel_error);	
end

% Get the spin-exchange cross-sections if they are there, but don't
% validate content yet (don't know the names of the other elements).
if ~exist('se', 'var')
	se = {element, NaN};		% Set its self-exchange to NaN at least.
end

serror = ['Spin exchange cross sections must be a cell array of ',...
	'two elements, the first being a cell array of (or a single) ', ...
	'element name(s), the second being a vector of (or a single) ', ...
	'spin exchange cross section(s).'];
if ~iscell(se) || length(se) ~= 2 || ~isnumeric(se{2}) || ...
		length(se{1}) ~= length(se{2})
	error(serror);
elseif iscell(se{1}) 
	if	any(cellfun(@(x)~ischar(x), se{1}))
		error(serror);
	end
elseif ~ischar(se{1})
	error(serror);		
end

% Everything is validated, dump it into a struct.
iso = struct('element', element, 'number', number, 'mass', mass, ...
	'spin', spin, 'qfs', qfs, 'abundance', abundance, ...
	'nd_kernel', nd_kern, 'QF', QF, 'se', se);



