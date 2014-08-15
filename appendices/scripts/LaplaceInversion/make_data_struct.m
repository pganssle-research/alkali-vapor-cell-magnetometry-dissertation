function ds = make_data_struct(x, y, z, std)
% Makes a data structure from either 1D or 2D data. If the third option is
% a 2D matrix, it is assumed that you're trying to make a 2D struct. If the
% third option is a scalar or missing, it is assumed that you are passing 
% the std.
% 
% 
% Inputs:
% 2D Mode:
% x:    First index vector (size n)
% y:    Second index vector (size m)
% z:    Matrix (size nxm)
% std: Standard deviation of data (Default: Estimated)
%
% 1D Mode:
% x:	Index vector  (size nx1)
% y:  Index vector	(size nx1)
% std:	Standard deviation of the data. (Default: Estimated)
%
% Outputs:
% ds -> Struct:
%       d.x -> x
%       d.y -> y
%       d.z -> z (2D Mode only)
%       d.std -> std
%
% Usage:
% ds = make_data_struct(x, y, z[, std]);
% ds = make_data_struct(x, y[, std]);

if(~isvector(x) || ~isvector(y))
   error('X and Y must be vectors.'); 
end

if(size(x, 1) < size(x, 2))
   x = x';
end

if(size(y, 1) < size(y, 2))
    y = y';
end

dims = 1;
if(exist('z', 'var'))
	sz = size(z);
	if(~isscalar(z))
		dims = 2;
		if(~(length(x) == sz(1)))
			 if(length(x) == sz(2))
				  z = z';
				  sz = size(z);
			 else
				 error('Z must be of size nxm'); 
			 end
		end

		if(~length(y) == sz(2))
			error('Z must be of size nxm'); 
		end	
	else
		std = z;
	end
end

if(~exist('std', 'var'))
	% Estimate the standard deviations. Comes from Eqn 49 of the referenced
	% paper in linvert/linvert2D.
	
	k = floor(length(y)/3);
	i = 0:(k-1);
	
	if(dims == 2)
		std1 = sum(sum((z(:, 3*i+1) - 2*z(:, 3*i+2) + z(:, 3*i+3)).^2));
		std1 = sqrt((1/(6*k*length(x)))*std1);

		k = floor(length(x)/3);
		i = 0:(k-1);

		std2 = sum(sum((z(3*i+1, :) - 2*z(3*i+2, :) + z(3*i+3, :)).^2));
		std2 = sqrt((1/(6*k*length(y)))*std2);

		std = std1*std2;
	else
		std = sum((y(3*i+1)-2*y(3*i+2) + y(3*i+3)).^2);
		std = (1/(6*k))*abs(std);
	end
end

if(dims == 2)
	ds = struct('x', x, 'y', y, 'z', z, 'std', std);
else
	ds = struct('x', x, 'y', y, 'std', std);
end


