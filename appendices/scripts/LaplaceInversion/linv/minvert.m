function [f, G, H] = minvert(c, Mc, K, alpha, id)
% Implements the unconstrained minimization from the referenced paper
% IEEE Transactions on Signal Processing, Vol. 50, No. 5 May 2002
%
% fr = max(0, K'c);

% This is the 'Kuhn-Tucker' condition.
% Here G(c) = K*min(0, Heavyside(diag(K'*c)))*K'
% The first two operators resolve to just K with the
% columns where K'*c <= 0 evalutes to true set to 0, so for speed
% we can use logical indexing and skip a matrix multiplication.
H = K;
H(:, K'*c <= 0) = 0;
H = H*K' + alpha*id;

G = H*c;	

f = 0.5*c'*G - c'*Mc;	% Eqn 31

if(nargout > 1)
	G = G - Mc;			% Eqn 32
end