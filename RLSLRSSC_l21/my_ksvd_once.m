function [ D ] = my_ksvd_once( Y, D, R ) %||Y-DR||_F^2  ||X-DW||
X = Y';
W = R';
[d n] = size(Y);
for k=1:n
    I = find(W(k,:));
    Ri = R(:,I) + D(:,k)*W(k,I);
    dk = Ri * W(k,I)';
    dk = dk/sqrt(dk'*dk);  % normalize
    D(:,k) = dk;
    W(k,I) = dk'*Ri;
    R(:,I) = Ri - D(:,k)*W(k,I);
end