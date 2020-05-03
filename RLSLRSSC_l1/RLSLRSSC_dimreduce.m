function [ P ] = my_dimreduce(P, Y, C, E, L, r, alpha1,alpha2)
% function [ P ] = my_dimreduce(P, Y,T, r, alpha1,alpha2)
% it solves the following problem
%-----------------------------------------
% min alpha1||PY-PYC-E||_F^2 + alpha2tr(PXLX'P')
% subject to: PP'=I;
%-----------------------------------------

% Delta = alpha1*tr((Y-PYC-P'*E)*(Y-PYC-P'*E)') + alpha2tr(PXLX'P');
Delta = alpha1*((Y-Y*C-P'*E)*(Y-Y*C-P'*E)') + alpha2*(Y*L*Y');
Delta = max(Delta,Delta');

[eigvec eigval] = eig(Delta);
eigval = diag(eigval);
[sort_eigval,sort_eigval_index]=sort(eigval,'ascend');
E=eigvec(:,sort_eigval_index);
P=E(:,1:r)';
end
