function [f, G, H] = minvert(c, Mc, K, alpha, id)
% Implements the unconstrained minimization from the referenced paper
% IEEE Transactions on Signal Processing, Vol. 50, No. 5 May 2002
%
% fr = max(0, K'c);

% This is the 'Kuhn-Tucker' condition.
% Here G(c) = K*max(0, Heavyside(diag(K'*c)))*K'
% The first two operators resolve to just K with the
% columns where K'*c <= 0 evalutes to true set to 0, so for speed
% we can use logical indexing and skip a matrix multiplication.

% Doing it this way, where you set H(:, K'*c) = 0, takes 0.34 ms/iteration,
% and K*diag(max(0, K'*c)>=0)*K' takes 0.68 ms/iteration.
Gc = K;
Gc(:, K'*c <= 0) = 0;	
Gc = Gc*K';

H = Gc + alpha*id;
G = H*c;	

f = c'*(0.5*G - Mc);	% Eqn 31

if(nargout > 1)
	G = G - Mc;			% Eqn 32
end