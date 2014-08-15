function y = diff_kernel(G, D, n, tau, gamma)
% Diffusion, for gradient G (G/cm) at diffusion coefficient D (1e-5 cm^2/s), 
% with n  CMPGs with interpulse spacing tau. Gamma is the gyromagnetic 
% ratio of the nucleus in Hz/Gauss.

if ~isrow(D)
	D = D';
end

if ~iscolumn(G)
	G = G';
end

A = -(2/3)*(2*pi*gamma)^2*n*tau.^3;
D = D*1e-5;
y = exp(A*(G.^2)*D);
