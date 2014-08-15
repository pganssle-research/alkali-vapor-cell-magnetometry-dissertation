function comps = add_lw_comp(name, mass, varargin)
% Gaseous component addition
% Last Updated: April 2013
% Author: Paul J. Ganssle,
%
% Function for adding new buffer gas components to the linewidth
% calculation. This generates a structure for each component which lw_opts
% can store for use by calculate_linewidth. The main properties of interest
% are how the gasses interact with one another and with the alkali metal
% atoms, so this generates a series of interaction matricies. They need not
% be complete for the program to operate, but be aware that this can
% introduce some errors.
%
% Of interest are the spin-destruction cross sections with each of the
% elements in lw_opts('get_element_names') [Default set: K, Rb, Cs], the
% diffusion coefficients of those elements in the gasses of interest at
% some reference temperature and pressure (D0, P0, T0).
%
% Adding an existing component merges the two records, favoring new
% records. Pass NaN for values you don't know if you need to, and pass a
% negative number to reset an erroneous value to NaN.
%
% Inputs:
% ------------------------------------------------------------------------
% name:			[String] The naem of the gas (e.g. N2, Ne, He)
% mass:			[Double] The mass of the gas in au (unified atomic mass)
% comps:        [Struct] Optional - a components structure to which this
%                  component will be added. If missing, a new structure
%                  will be generated.
%
% The other inputs are important because they begin to fill out a grid of
% information. Elements are not checked for existence until the struct is
% passed to lw_opts. The properties are structured as:
% {'Element', {'prop1', val1, 'prop2', val2, ...} where 'Element' is the
% name of an element, 'prop1' to 'propN' is the name of a valid property,
% and val can be any valid value - generally either a scalar or vector.
% The valid properties are:
%
% 'sd':			'sd', [scalar] - Sets the spin destruction cross-section
%					 between this component and 'Element'.
%
% 'se':			'se', [scalar] - Sets the spin exchange cross-section
%					 between this component and 'Element'.
%
% 'sq':			'sq', [scalar/vector] - Sets the spin quenching cross-
%                   section. If it's a vector, it should be a 2-vector
%                   [P1/2, P3/2].
%
% 'D0':			'D0', [D0, P0, T0] - Sets the reference point for
%					Element's diffusion coefficient (D0) at a reference
%					pressure (P0) and temperature (T0). If a scalar is
%					passed to 'D0' and D0 has not been set before, P0 is
%					set to 1 atm and T0 is set to 298.15K. If P0 and T0
%                   have not been set before, only D0 is replaced.
%
% ------------------------------------------------------------------------
%
% Output:
% ------------------------------------------------------------------------
% comps:		A struct containing all the information. If 'comps' has
%               been passed as an argument, the component is added into
%               the components structure.
% ------------------------------------------------------------------------
%
% Prototypes:
% comps = add_lw_comp(number, mass, props);
% comps = add_lw_comp(number, mass, comps, props);

% Input validation
if ~ischar(name)
    error('Gas name must be a string.');
end

if ~isnumeric(mass)
    error('Mass must be a numeric value in atomic mass units.');
end

% Check if a previous version has been loaded.
if(isstruct(varargin{1}))
    comps = varargin{1};
    varargin(1) = [];
else
    comps = struct();
end

% Read in previously loaded values, if available.
if isfield(comps, name)
    s = comps.(name);
end

% Make a new structure if one didn't exist.
s.mass = mass;

% Check for the 4 entries.
for ii = 1:length(varargin)
    entry = varargin{ii};
    if ~iscell(entry)
        continue;
    end
    
    if ischar(entry{1})
        elem_name = entry{1};
    else
        warning('Unnamed element passed to arguments. Ignoring.');
        continue;
    end
    
    if length(entry) ~= 2
        warning(['Invalid argument passed to ', elem_name, ' - value', ...
            ' is too long or too short. Ignoring.']);
        continue;
    end
    
    if ~iscell(entry{2})
        warning(['Invalid argument passed to ', elem_name, ' - value', ...
            ' must be a cell of type {''prop'', value, ...}. Ignoring']);
        continue;
    end

    vals = entry{2};
    valid_props = {'sd', 'sq', 'se', 'D0'};
    if ~any(cellfun(@(x)any(strcmpi(x, vals(1:2:end))), valid_props))
        warning(['No valid properties passed to ', elem_name, ...
            ', ignoring.']);
    end
    
    if isfield(s, elem_name)
       elem = s.(elem_name); 
    end
    
    % sd
    sd = parse_varargin('sd', vals); % Don't check types yet.
    if ~isempty(sd)
        if ~isnan(sd) && (~isnumeric(sd) || ~isscalar(sd))
            warning(['Invalid spin destruction cross-section was passed ', ...
                ' for element ', elem_name, '. Ignoring value.']);
        elseif ~isnan(sd) || ~isfield(elem, 'sd')
            elem.sd = sd;
        end
    end
    
    % se
    se = parse_varargin('se', vals); % Don't check types yet.
    if ~isempty(se)
        if ~isnan(se) && (~isnumeric(se) || ~isscalar(se))
            warning(['Invalid spin exchange cross-section was passed ', ...
                ' for element ', elem_name, '. Ignoring value.']);
        elseif ~isnan(se) || ~isfield(elem, 'se')
            elem.se = se;
        end
    end
    
    % sq
    sq = parse_varargin('sq', vals); % Don't check types yet.
    if ~isempty(sq)
        if ~any(isnan(sq)) && (~isnumeric(sq) || length(sq) > 2)
            warning(['Invalid spin quenching cross-section was ', ...
                'passed for element ', elem_name, '. Ignoring value.']);
        elseif ~isfield(elem, 'sq')
            elem.sq(~isnan(sq)) = sq(~isnan(sq));
        end
    end
    
    % D0, P0, T0
    D0 = parse_varargin('D0', vals); % Don't check types yet.
    if ~isempty(D0)
        if ~isnan(D0) && (~isnumeric(D0) || (~isscalar(D0) ...
                && length(D0) ~= 3))
            warning(['Invalid diffusion coefficient was passed ', ...
                ' for element ', elem_name, '. Ignoring value.']);
        elseif ~isnan(D0) || ~isfield(s, 'D0')
            if isscalar(D0)         % Should be a 3-vector, or use STP.
                P0 = 760;           % 1 atm in torr 
                T0 = 298.15;        % RT in K.
            else
                P0 = D0(2);
                T0 = D0(3);
                D0 = D0(1);
            end
            
            elem.D0 = D0;   elem.P0 = P0;   elem.T0 = T0;
        end
    end
    
    if ~isempty(elem)
       s.(elem_name) = elem; 
    end
end

comps.(name) = s;