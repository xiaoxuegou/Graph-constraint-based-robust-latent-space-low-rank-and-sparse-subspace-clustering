function [ A ] = RLSLRSSC_initial_P(X, L, r)

GL = X*L*X';
GL = max(GL,GL');
[eigvec, eigval] = eig(GL);

eigval = diag(eigval);%get the main diagonal element
[~,sort_eigval_index] = sort(eigval,'ascend');%sort the eign value orderly
T = eigvec(:,sort_eigval_index);%sort the eign vector according to eign value
A = T(:,1:r)';
end


