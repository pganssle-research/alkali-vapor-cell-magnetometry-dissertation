function val = parse_varargin(keys, varg, types, bound)
% Variable Argument Parser
% Last Updated: April 2013
% Author: Paul J. Ganssle, Pines Lab, UC Berkeley, Berkeley, CA
%
% Function for getting the values from a key in variable input arguments,
% as in many Matlab functions, with the form, 'Property', 'Value'. This
% searches a cell array for a specific property (key), or a cell array of
% keys. Optional arguments allow for the type and value to be constrained.
% Returns only the first matching argument for each
%
% Inputs:
% ------------------------------------------------------------------------
% keys:     Can be either a string with the key or a cell containing all
%			variations on the name.
%
% varg:		The varargin from the calling function. Should be a cell.
%
% types:	Allowed types for the given key. Default is -1 (not specified).
%			Otherwise pass a cellstr with any of the following:
%					'num'
%					'vector'
%					'matrix'
%					'cell'
%					'str'
%					'bool'
%			For multiple values of 'types', first the result is tested as a
%           num, vector, matrix, cell, str, then a bool. If a cell array is
%           passed to 'keys', use a cell array of cells to specify
%           individual type bounds, (must have 1 entry per entry in keys),
%           otherwise this value is applied to all entries.
%
% bound:	Bounds for the values. Only really valuable for nums, has the
%			form [lb, ub] Default = -1 (no bounds). For a cell array input
%			to keys, to use different values for each parameter, you must
%			provide a cell array of vectors, otherwise the same values will
%			be applied to each element in keys.
% ------------------------------------------------------------------------
%
% Outputs:
% ------------------------------------------------------------------------
% val:      Assumiung no errors, returns the value associated with the
%           property. If the property is not present, returns an empty
%           vector.
% ------------------------------------------------------------------------
%
% val = parse_varargin(keys, varg, types, bound);

% Useful for later - anonymous functions for converting strings
% describing boolean values.
istruestr = @(lstr)any(strcmpi(lstr, {'t', 'true', 'y', 'yes'}));
isfalsestr = @(lstr)any(strcmpi(lstr, {'f','false', 'n', 'no'}));

% Default values for whether or not each type check is necessary.
visnum = false;     visvec = false;     vismat = false;
viscell = false;    visstr = false;     visbool = false;

hastype = false;     hasbound = false;

if exist('types', 'var') && (iscell(types) || ischar(types))
    if(~iscell(types))
        types = {types};
    end
    
    types = types(cellfun(@(x)any(strcmp({'num', 'vector', 'matrix', 'cell', 'str', 'bool'}, x)), types));
    if ~isempty(types)
        hastype = true;
        visnum = any(strcmpi('num', types));
        visvec = any(strcmpi('vector', types));
        vismat = any(strcmpi('matrix', types));
        viscell = any(strcmpi('cell', types));
        visstr = any(strcmpi('str', types));
        visbool = any(strcmpi('bool', types));
    end
end

if exist('bound', 'var') && isvector(bound) && length(bound) == 2
    hasbound = true;
end

if ~iscell(keys)
    keys = {keys};
end

% Check all the keys against the varargin
vals = zeros(size(varg));
vals(1:2:end) = cellfun(@(x)any(strcmp(x, keys)), varg(1:2:end));
ints = find(vals);

if isempty(ints)
    val = [];
    return
end

% Check for the first valid value.
for i = ints
    val = [];
    out = varg{i+1};
    
    if isnumeric(out) && hasbound
        if any(out < bound(1)) || any(out > bound(2))
            warning('Out of bounds value found.')
            continue
        end
    end
    
    if hastype
        if (visnum && isscalar(out)) || ...
                (visvec && isvector(out)) || ...
                (vismat && ismatrix(out)) || ...
                (viscell && iscell(out)) || ...
                (visstr && ischar(out))
            val = out;
        elseif visbool
            % Even if it's a char, it might be a logical string.
            if ischar(out)
                truestr = istruestr(out);
                
                if truestr || isfalsestr(out)
                    val = truestr;
                end
            end
            
            if isnumeric(out)
                val = logical(out);
            end
        end
    else
        val = out;
    end
    
    % Found the first valid one.
    if ~isempty(val)
        return
    end
end