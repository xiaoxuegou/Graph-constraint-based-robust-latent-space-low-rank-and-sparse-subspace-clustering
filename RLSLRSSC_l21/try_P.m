function [history] = try_P()
% it solves the following problem
% min 0.5||PY-PAC||^2_F + lambda1||C||_* + lambda2||C||_1 + lambda3||Y-P'PY||^2_F
% subject to: C^T*1_n=1_n
%-----------------------------------------
%P: D*r   Y:D*N   A:D*k   C:k*N
% ITERMAX = 3;


Y = rand(9,20);
ADAPTIVE_ON = 1;
[d n] = size(Y);
 r = 3; 
%  k = ceil(0.5*n);  
 k = 5;  
 lambda1 = 1;
 lambda2 = 1;
 lambda3 = 1;

%%%%%Initialize:
%randID = randperm(n);
% randID = 1:n;
% D0 = dictnormalize(Y(:, randID(1:k)));
% D0 = dictnormalize(rand(d,k));
D0 = rand(d,k);
C0 = my_ALM_noisyLRSSC( D0, Y, lambda1, lambda2, ADAPTIVE_ON); 
P0 = 0;

A = D0;
C = C0;
P = P0;

W = pinv(Y) * A *  C;
P = my_dimreduce(Y, W, r, 0.5, lambda3);%step1
% PA =  compute_dictionary_ODL( P*Y, P*A, C );%step2
% 
% C = my_ALM_noisyLRSSC( PA, P*Y, lambda1, lambda2, ADAPTIVE_ON);%step3
% A = pinv(P)*PA;
dual_trans = sum(sum((Y-P'*P*Y).^2));
dual_pinv = sum(sum((Y-pinv(P)*P*Y).^2));

history.Y = Y;
history.W = W;
history.C = C;
return;