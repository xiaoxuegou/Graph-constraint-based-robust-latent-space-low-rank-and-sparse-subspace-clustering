function [D] = UpdateD_ALM(X,D,C,lambda,eta)
[d,K] = size(D);
D0 = D;
D1 = D;

u = eta/5;

Lambda = zeros(d,K);
max_u = 1e30;
rho0 = 1.1;
epsilon = 1e-8;
max_iter = 300;
iter = 1;
converge = 0;
while ~converge
    %update D1
    D1_new = sigma_soft_thresh(Lambda/u+D0,eta/u);
    D1_new = D1_new*diag(1./sqrt(diag(D1_new'*D1_new)));
    
    %update D0
    D0_new = (2*lambda*X*C'+u*D1_new-Lambda)*pinv(lambda*(C*C')+u*eye(K));
    D0_new = D0_new*diag(1./sqrt(diag(D0_new'*D0_new)));
    
    % D0_new = D1_new;
    % A = lambda*(C*C') + (u/2)*eye(K);
    % B = 2 * lambda * X * C' + u * D0_new - Lambda;
    %can we do not have convergence analysis ?
    % temp = D0_new;
    % c = 0;
    % for i = 1:K
       %  ui = (B(:,i) - D0_new * A(:,i)) / A(i,i) + D0_new(:,i);
       % D0_new(:,i) = ui / max(norm(ui),1);
    % end
    
    res_obj = norm(X-D0_new*C,inf);
    dual_res = norm(D0_new - D1_new,inf);
    
    D0 = D0_new;
    D1 = D1_new;
    
    if dual_res <= epsilon || iter > max_iter
        converge = 1;
    end
    if mod(iter,50) == 0
         fprintf('%d th iteration compelte, res_obj= %f, dual_res=%f\n',iter,res_obj,dual_res);
    end
    
    Lambda = Lambda + u*(D0-D1);
    u = min(rho0*u,max_u);
    iter = iter + 1;
end
D = D0;
end