function [D,C,E] = my_inexact_ALM(PX,D,lambda1,lambda3)
%min ||C||_*+lambda1||C||_1+lambda3||E||_l
[d,N] = size(PX);
k = size(D,2);
J = zeros(k,N);
C1 = zeros(k,N);
C2 = zeros(k,N);
E = zeros(k,N);

eta = 1e-8;
u = 1e-6;
max_u = 1e6;
rho = 1.1;

Lambda1 = zeros(k,N);
Lambda2 = zeros(N,1);
Lambda3 = zeros(k,N);
Lambda4 = zeros(k,N);

converge = 0;
%%
while ~converge
    %update J
    J = pinv(u*D'*D+2*u*eye(k)) * (D'*Lambda1-ones(N,1)'*Lambda2-Lambda3-Lambda4+u*D'*PX-u*D'*E+u*ones(N,N)+u*(C1+C2));
    %%
    %update C1
    C1 = sigma_soft_thresh(J+Lambda3/u,1/u);

    %%
    %update C2
    C2 = soft_thresh(J+Lambda4/u,lambda1/u);
    %%
    %update D
    D = (Lambda1/u+PX-E)*pinv(J);
    %%
    %update E
    E = (Lambda1+u*PX-u*D*J)/(2*lambda3+u);
    %%update Lambda
    Lambda1 = Lambda1+u*(PX-D*J-E);
    Lambda2 = Lambda2+u*(J*ones(N,1)-ones(N,1));
    Lambda3 = Lambda3+u*(J-C1);
    Lambda4 = Lambda4+u*(J-C2);
    %%
    %update u
    u = min(rho*u,max_u);
    %%
    %judge convergence
    pre_1 = norm(PX-D*J-E,'inf');
    pre_2 = norm(J*ones(N,1)-ones(N,1),'inf');
    pre_3 = norm(J-C1,'inf');
    pre_4 = norm(J-C2,'inf');
    while min(pre_1,pre_2,pre_3,pre_4) < eta
        converge = 1;
        C = C2;
    end
end
end