function out = process_data(in, poly_order, cut, apod_lb)
% Does the most basic data processing - polynomial subtraction, truncation
% and zero-packing
%
% out = process_data(in[, poly_order, cut, apod_lb])
if(~exist('poly_order', 'var'))
    poly_order = -1;
end

if(~exist('cut', 'var'))
    cut = -1;
end

out = in;

% In case you want to make adjustments.
if(~isfield(out, 'odata'))
    out.odata = out.mdata; % Save the old data
    out.mdata = sub_poly(out.mdata, poly_order, 0, cut+1);
else
    if(cut >= 0 || poly_order >= 0)    
        out.mdata = sub_poly(out.odata, poly_order, 0, cut+1);
    end
end

% Cut and zero pack
if(isscalar(cut) && cut > 1)
    out.mdata = circshift(out.mdata, -cut);
    out.mdata((end-cut+1):end, :) = 0;
end

% Apodize
if(exist('apod_lb', 'var') && apod_lb > 0)
    out.mdata = apodize(out.mdata, out.t, apod_lb);
end