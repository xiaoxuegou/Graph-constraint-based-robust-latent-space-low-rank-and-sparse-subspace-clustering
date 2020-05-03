function [L] = initial_laplace(X,options)
if nargin < 2
    options = [];
end
K = 8;
if isfield(options,'K')
    K = options.K+1; %different smaples corresponding to  differnet sigma 
end
k = 4;%4-10,according to the size of dataset
if isfield(options,'k')
    %the parameter of knn
    k = options.k; 
end

%% initial the sigma matrix
N = size(X,2);
single = zeros(N);
temp_sigma = zeros(N,1);
for i = 1:N
    for j = 1:N
        temp_sigma(j,1) = norm(X(:,i)-X(:,j),2);
    end
    temp_sigma = sort(temp_sigma,'ascend');
    single(i) = temp_sigma(K);
%         single(i) = mean(temp_sigma);
end
single_matrix = single * single';
    
%% compute the unnormlized laplace matrix
W = zeros(N,N);
C = zeros(N,N);
for i = 1:N
    for j = i:N
        dif = X(:,i)-X(:,j);
        W(i,j) = exp(-sum(dif.*dif) / (single_matrix(i,j)));
    end
end
W = W + W';
W = W - diag(diag(W));
for i = 1:N
    [~,index] = sort(W(i,:),'descend');
    C(i,index(1:k)) = W(i,index(1:k));
end
%compute L
C = max(C,C');
D = sum(C,2);
D = diag(D);
L = D - C;
end