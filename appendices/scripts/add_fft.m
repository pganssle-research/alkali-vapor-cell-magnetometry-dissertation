function out = add_fft(in, poly_order, cut, apod_lb)
% Give this a structure like those generated from the pulse controller,
% this will return a structure with ffts added to it.

if(~exist('poly_order', 'var'))
    poly_order = -1;
end

if(~exist('cut', 'var'))
    cut = -1;
end

if(~exist('apod_lb', 'var'))
    apod_lb = -0.5;
end

% Do the basic data processing.
if(poly_order >= 0 || cut > 0 && apod_lb > 0)
	out = process_data(in, poly_order, cut, apod_lb);
else
	out = in;
end

% Apply the fourier transform.
sr = out.prog.sr;
np_fft = 2^(ceil(log2(size(out.mdata, 1)))+1);

f = linspace(0, sr/2, np_fft/2);
s = fft(out.mdata(:, :), np_fft);
s = 2*s/size(out.mdata, 1);

sd = num2cell(size(out.mdata));

out.f = f;
out.fft = s(1:(np_fft/2), :);

if(length(sd) > 1)
    out.fft = reshape(out.fft, np_fft/2, sd{2:end});
end




