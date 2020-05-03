function [D] = metaface(P,Y,A,C)
X = P*Y;
[K] = size(A,2);
[d] = size(X,1);
Z = zeros(d,K);
for i = 1:K
    Z = Z + A(:,i)*C(i,:);
end
for i = 1:K
    Z = Z - A(:,i)*C(i,:);
    A(:,i) = (Z * C(i,:)') / norm(Z * C(i,:)');
    Z = Z + A(:,i)*C(i,:);
end
D = A;
end