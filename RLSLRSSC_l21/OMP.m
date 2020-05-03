function [S] = OMP(X,D,L)
% X: input signal
% D: dictionary, each column is an atom. Note: atoms much be normalized.
% L: number of non-zero entries in output
%   output, a matrix that stores the coefficient of each atom
S = [];
n = size(X,2);
for i = 1:n
    t = omp(X(:,i),D,L);
    S = [S t];
end

end

function [t] = omp(x, D, L)
% x: input signal
% D: dictionary, each column is an atom. Note: atoms much be normalized.
% L: number of non-zero entries in output
%   output, a vector that stores the coefficient of each atom

dim = size(D,2);
t = zeros(dim,1);
r = x; % residual
selected_atom = [];
for iter = 1 : L
    dot_p = D' * r;
    [~, i] = max(abs(dot_p));
    selected_atom = [selected_atom i];
    temp = D(:, selected_atom);
    a = pinv(temp) * x;
    r = x - temp * a;
end
t(selected_atom) = a;
end