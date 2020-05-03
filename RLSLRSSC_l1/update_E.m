function [E] = update_E(P,Y,C,lambda,lambda4,E_norm)
X = P*Y-P*Y*C;
thresh = lambda4/(2*lambda);
if E_norm == 1
    E = soft_thresh(X,thresh);
else
    %l21_norm
    [d,N] = size(X);
    E = zeros(d,N);
    for i = 1:size(X,2)
        i_norm = norm(X(:,i));
        if i_norm > thresh
            E(:,i) = X(:,i)*(i_norm-thresh)/i_norm;
        end
    end
end
end