function [ C1, history ] = ADMM_RLSLRSSC( A, B,lambda, lambda1,lambda2, ADAPTIVE_ON)%A,B
%ALM_LowRankSSC Summary of this function goes here
%   Detailed explanation goes here
% sigma is the estimated noise level
% lambda defines the trade-off between nuclear norm and one-norm

% it solves the following problem
% min lambda||B-AC||^2_F + lambda1||C||_* + lambda2||C||_1
% subject to: C^T*1_k=1_n,diag(C)=0

%A:d*k(initial dictionary)    B:d*n(data samples)    C:k*n    I_k:k*1    I_n:n*1    Lambda1:1*n    Lambda2:k*n    Lambda3:k*n
%-----------------------------------------
% lambda1=2*sigma*beta/n
% lambda2=2*sigma*(1-beta)/n

lambda = lambda * 2;
[d,n]=size(B);%d*n
[trash,k]=size(A);%d*k

% t=1+lambda;
% beta=1/t;
% beta2=1-beta;
% 
% lambda1=beta/sqrt(2*log(n))*sigma;
% lambda2=beta2/sqrt(2*log(n))*sigma;



mu=0.2;
mu1=mu;
mu2=mu;
mu3=mu;
% the ratio between mu1, mu2, mu3 can be tuned
% but once these are changed once, one matrix inversion need to be computed


% optimal variables
rho0 = 1.1;% this one can be tuned
rho=1;

epsilon=0.1;
eta=epsilon*5*n;
max_mu=1000;

tol_abs=1e-5;
tol_rel=1e-3;

max_iter=400;
iter=0;

%initialze
C1=zeros(k,n);
C2=zeros(k,n);
J=zeros(k,n);
I_k=ones(k,1);
I_n=ones(n,1);
Lambda1=zeros(1,n);
Lambda2=zeros(k,n);
Lambda3=zeros(k,n);
fun_val=0;


%caching a few quantities
ATA=(A'*A); %
ATB=(A'*B);
invMat=pinv(lambda*ATA+mu1*ones(k,k)+(mu2+mu3)*eye(k));% k*k

% diag_idx=logical(eye(n));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fix me, delete?
converge=0;
while ~converge
    if rho~=1
        invMat=pinv(lambda*ATA+mu1*ones(k,k)+(mu2+mu3)*eye(k));
    end
    

    J_new=invMat*(lambda*ATB+mu1*ones(k,n)+mu2*C2+mu3*C1-I_k*Lambda1-Lambda2-Lambda3);%
    C2_new=soft_thresh(J_new+Lambda2/mu2,lambda2/mu2);
%     C2_new(diag_idx)=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fix me, delete?
    temp = J_new+Lambda3/mu3;
    C1_new=sigma_soft_thresh(temp,lambda1/mu3);
    
    
    %primal and dual residuals
    d_res=mu2*norm(C1_new(:)-C1(:))+mu3*norm(C2_new(:)-C2(:));
    %pDist1=J_new'*I_n-I_n;
    pDist1=I_k'*J_new-I_n';%J_new'*I_k-I_n
    pDist2=J_new-C2_new;
    pDist3=J_new-C1_new;
    p_res=norm(pDist1(:))+norm(pDist2(:))+norm(pDist3(:));
    
    %update
    C1=C1_new;
    C2=C2_new;
    J=J_new;
    Lambda1=Lambda1+mu1*pDist1;
    Lambda2=Lambda2+mu2*pDist2;
    Lambda3=Lambda3+mu3*pDist3;
    
    if ADAPTIVE_ON
        rho=min(max_mu/mu2,rho0);
        if d_res*eta/norm(B(:))>epsilon
            rho=1;
        else
            mu1=rho*mu1;
            mu2=rho*mu2;
            mu3=rho*mu3;
        end
    end
    
    
   % It is pointless to know function value
   % cur_fun_val=lambda1*sum(svd(C1))+lambda2*norm(C1(:),1);
   % res_func=abs(cur_fun_val-fun_val); 
     
   iter=iter+1;
   if mod(iter,20)==0
    fprintf('%d th iteration compelte, P_res= %f, D_res=%f\n',iter,p_res,d_res);
%     fprintf('    Feasibility is :[%f, %f] mu=[%.2f, %.2f]\n',norm(pDist2(:)),norm(pDist3(:)),...
%         mu2,mu3);
   end
   history(iter)=p_res;
   
   
   if( p_res<n*tol_abs+tol_rel*max(norm(C1,'fro'),norm(J,'fro')) && d_res <n*tol_abs+tol_rel*norm([Lambda1(:);Lambda2(:);Lambda3(:)]))
%   if( p_res<n*tol_abs && d_res <n*tol_abs)
       converge=1;
       fprintf('Converged.\n')
   end
   if iter>=max_iter
   converge=1;
   fprintf('Maximum iteration reached. Quit.\n')
   end

   
   %fun_val=cur_fun_val;
end
end