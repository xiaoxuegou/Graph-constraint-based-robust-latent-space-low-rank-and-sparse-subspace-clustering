function [fea] = nomalize_each_vector_to_unit(data)
fea = data;
%Nomalize_each_vector_to_unit
%===========================================
[nSmp,nFea] = size(fea);
for i = 1:nSmp
     fea(i,:) = fea(i,:) ./ max(1e-12,norm(fea(i,:)));
end
%===========================================
return;