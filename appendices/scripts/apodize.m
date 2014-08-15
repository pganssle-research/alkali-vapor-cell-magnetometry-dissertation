function out = apodize(data, t, lb)
% Function that generates apodized data in 1 dimension
% data = data you want in
% t = time vector in seconds
% lb = line broadening in hz
%
% size(data, 1) must be the same size as t.
% Applies along only the first dimension.
%
% Usage: out = apodize(data, t, lb)

if(size(data, 1) ~= length(t))
   error('Incommensurate size between data and time vectors.'); 
end

if(size(t, 1) ~= length(t))
   t = t'; 
end

% Generate the exponentials
e_t = exp(-t/lb);

% Apply it to the data
data(:, :) = cell2mat(arrayfun(@(x)e_t.*data(:, x), ...
                      1:length(data(1, :)), ...
                      'UniformOutput', false));

out = data;