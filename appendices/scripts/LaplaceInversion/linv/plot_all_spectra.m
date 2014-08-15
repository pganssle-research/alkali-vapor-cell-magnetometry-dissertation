function plot_all_spectra(s, n, start, s1, s2)
% Function for plotting all the spectra output by linvert or linvert2D, in
% case you want to see the history of all the alphas and such.
%
% Inputs;
% s - The structure output from linvert or linvert2D
% n - The number of alphas you want to see. Pass <= 0 for default [all]
% start - The alpha index you want to start at. Pass <= 0 for default [all]
% s1 - Subplot dimension 1 (default = ceil(sqrt(n))
% s2 - Subplot dimension 2 (default = ceil(n/s1))
%
% Usage:
% plot_all_spectra(s, n, start, s1, s2);

if(~exist('n', 'var') || n <= 0 || n >= length(s.f))
	n = length(s.f);
end

if(~exist('start', 'var') || start <= 0 || start+n > length(s.f)+1)
	start = 1;
end

if(~exist('s1', 'var'))
	s1 = ceil(sqrt(n));
end

if(~exist('s2', 'var'))
	s2 = ceil(n/s1);
end

n = max(n, s1*s2);
if(isvector(s.f{1}))
	dims = 1;
else
	dims = 2;
end

for i = 0:(n-1)
	subplot(s1, s2, i+1);

	if(dims == 2)
		contour(s.t1, s.t2, s.f{i+start}');
	else
		plot(s.t, s.f{i+start});
	end
	
	title(sprintf('%d: \\alpha = %3.2f', i+start, s.alpha(i+start)));
end