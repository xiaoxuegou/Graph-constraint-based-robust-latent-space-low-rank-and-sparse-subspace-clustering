function [ P ] = dim_reduce(P, Y,T, r, alpha1,alpha2,C)
% function [ P ] = my_dimreduce(P, Y,T, r, alpha1,alpha2)
% it solves the following problem
%-----------------------------------------
% min alpha1||PY-T||_F^2 + alpha2||Y-P'PY||_F^2
% subject to: PP'=I;
%-----------------------------------------

% Delta = alpha1*((Y-P'*T)*(Y-P'*T)') - alpha2*(Y*Y');
Delta = alpha1*(Y-Y*C-P'*T)*(Y-Y*C-P'*T)' - alpha2*(Y*Y');

Delta(isnan(Delta)) = 0;
Delta(isinf(Delta)) = 0;
Delta = (Delta + Delta') / 2;

[eigvec eigval]=eig(Delta);
eigval=sum(eigval,1);
[sort_eigval,sort_eigval_index]=sort(eigval,'ascend');
T=eigvec(:,sort_eigval_index);
P=T(:,1:r)';
end
