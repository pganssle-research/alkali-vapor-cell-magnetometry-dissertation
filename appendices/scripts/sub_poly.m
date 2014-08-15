function out = sub_poly(in, order, windowed, cut)
% Function that subtracts a polynomial offset from your data.
%
% Feed this a 1D or 2D matrix, tell it how many points you want to cut off, and
% it will fit those points to a polynomial and apply said polynomial to the
% output.
%
% If windowed evaluates to logical-true, then the first two dimensions will
% be reshaped before and after evaluation.
%
% Usage: out = sub_poly(in, order, windowed, cut);

if(~exist('windowed', 'var'))
    windowed = 0;
end

s_orig = size(in);
s_orig_cell = num2cell(s_orig);

if(windowed)
    if(length(s_orig) > 2)
        in = reshape(in, prod(s_orig(1:2)), s_orig_cell{3:end});
    elseif(length(s_orig)>1)
        in = reshape(in, prod(s_orig(1:2)), 1);
    end
end

s = size(in);

if(~exist('order', 'var') || order < 0)
	order = 1;
end

if(~exist('cut', 'var') || (isscalar(cut) && cut < 1) || isempty(cut))
    cut = 1;
end

if order > 99
	order = 99;
end

points = (1:s(1))'; % Initialize the "t" vector to have all the points

% Get the polynomial fits
if(isscalar(cut))
    p = cell2mat(arrayfun(@(x)polyfit(points(cut:end), ...
                                      in(cut:end, x), ...,
                                      order), ...
                          1:prod(s(2:end)), ...
                 'UniformOutput', false)');
else
    p = cell2mat(arrayfun(@(x)polyfit(points(cut), ...
                                      in(cut, x), ...
                                      order), ...
                          1:prod(s(2:end)), ...
                'UniformOutput', false)');
end

% Get the output array
out = cell2mat(arrayfun(@(x)in(:, x)-polyval(p(x, :), points), ...
                        1:prod(s(2:end)), ...
               'UniformOutput', false));

% Needs to be reshaped
if(windowed)
    if(length(s_orig) > 1)
        out = reshape(out, s_orig_cell{:});
    end
else
    s = num2cell(s);
    out = reshape(out, s{:});
end