function y = exp_decay(x, tau1)
% Simple exponential fit.

y = cell2mat(arrayfun(@(t)exp(-x/t), tau1, 'UniformOutput', false));