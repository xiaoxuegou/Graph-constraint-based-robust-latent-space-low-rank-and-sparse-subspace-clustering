function [D] = DCT_D(Y,K)
%DCT dictionary
Pn=ceil(sqrt(K));
bb=ceil(sqrt(size(Y,1)));
DCT=zeros(bb,Pn);
for k=0:1:Pn-1,
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0
        V=V-mean(V);
    end
    DCT(:,k+1)=V/norm(V);
end
DCT=kron(DCT,DCT);
D = DCT(1:size(Y,1),1:K);
end