function [dydx, x, dydx2, x2] = derivative(x, y, n, dim, inter)
% Take the nth derivative of a given function x with respect to the
% function y along dimension dim. Passing logical true to the argument
% 'inter' returns a spline-interpolation of the derivative function to the
% original data. This is only available for 1D data sets at the moment.
%
% Passing empty set to x gives a size-preserving version of diff(y, n, dim)
%
% Default n = 1;
% Default dim = 1 (or for vectors the primary axis);
% Default inter = true for 1D, false for ND;
%
% Outputs:
% dydx:		The nth derivative of y with respect to x, spline-interpolated
%				if available and 'inter' is set to true.
%
% x:   		The x function with size matched to dydx. For interpolated 
%				outputs, this is the same as the input x.
%
% dydxns:	The dydx function without spline interpolation. For inter =
%				false, this is the same as dydx.
%
% xns:		The x function without spline interpolation. For inter = false,
%				this is the same as dydx.
%
% [dydx, x, dydxns, xns] = derivative(x, y, [n, dim, inter])

% These must be defiend before we check for optional inputs.

dydx = y; % 0th derivative.

% What dimension should we vary along
if ~exist('dim', 'var') || isempty(dim) || dim < 0 || dim > length(size(y))
	dim = 1;
end

nd = length(size(dydx));
if nd == 2
	data1d = isvector(dydx);
end

if data1d
	if isrow(dydx)
		dim = 2;
	else
		dim = 1;
	end
end

if isempty(x)
	x = 1:size(y, dim);
elseif length(x) ~= size(y, dim)	
	error('x and y must have the same size.');
end
xx = x;

if ~exist('n', 'var') || n <= 0
	n = 1;
elseif n == 0
	dydx2 = dydx;
	x2 = xx;
	return;
end

nl = ones(1, length(size(y)));
nl(dim) = length(xx);
xx = reshape(xx, nl);

% Whether to spine-interpolate the output.
if ~exist('inter', 'var')
	inter = true;
else
	try
		inter = logical(inter);
	catch %#ok
		warning('Value for ''inter'' cannot be converted to logical. Using default of true');
		inter = true;
	end
end	

% Used for the eval in choosing x2 - don't know of a better way to do this
% than eval.
% nd = length(size(xx));
% pstr = ['xx(' repmat(':,', 1, dim-1) '%s,' repmat(':,', 1, nd-dim)];
% pstr = [pstr(1:(end-1)), ')'];
% pstr = ['(' sprintf(pstr, '2:end') '+' sprintf(pstr, '1:(end-1)') ')/2;'];

% Here's one:
% lcell = arrayfun(@(x)1:x, size(xx));
% lc1 = lcell;
% lc2 = lcell;
% lc1{dim} = 1:(size(xx, dim)-1);
% lc2{dim} = 2:(size(xx, dim));

% Doesn't matter anyway, this only ever takes the derivative along a single
% dimension, no need to fuck with all that.
rs = size(y);
rs(dim) = 1;

% Recursively calculate the derivative. Rather than using the recursive
% functionality of diff(), the recursion is done here so that the diff(x)
% step can be added in each time.
for i = 1:n
	dy = diff(dydx, 1, dim);
	dx = diff(xx, 1, dim);
	
% 	xx = eval(pstr);
% 	xx = (xx(lc1) + xx(lc2))/2;
	
	xx = (xx(1:(end-1))+xx(2:end))/2;
	
	dydx = dy./repmat(dx, rs);
end

dydx2 = dydx;
x2 = xx;

% Spline interpolation should be done at the end, as it introduces some
% model error, and that model error may be largest in the low-order 
% derivatives (i.e. spline extrapolation from x^3 will be more error-prone
% than spline extraplation from d^3(x^3)/d(x)^3 = 6.
if inter
	permin = 1:length(size(dydx));
	permin([2, dim]) = permin([dim, 2]);
	dydx = ipermute(spline(xx, permute(dydx, permin), x), permin);
else
	x = xx;
end