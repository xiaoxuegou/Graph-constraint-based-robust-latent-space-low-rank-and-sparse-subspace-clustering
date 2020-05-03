function [D,X] = Initializa_dictionary(Y,K,iter_max,S)
d = size(Y,1);
F = dictnormalize(DCT_D(Y,K)); %initial the F, and p=K?
if(nargin < 4)
   %S = ceil(size(Y,1)/20); 
   S=20;
end
if(nargin < 3)
   iter_max = 50; 
end
X = OMP(Y,F,S);
D = F; 
count = 0;
while count <= iter_max
    Delta = Y * X' * D' + eps;
    [U,~,V] = svd(Delta,'econ');
    %Q = V * U';
    Q = U * V';
    D = Q * D;
    X = OMP(Y,D,S);
    if(mod(count,10) == 0)
        fprintf('Iter %d, residual: %f\n', count, norm(Y-D*X,'fro'));
    end
    count = count + 1;
end
end