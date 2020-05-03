function [pca_vector,data] = myPCA(X, final_dimension, times)
%get the pca projection vector
if (~exist('pcaRatio','var'))
    pcaRatio = 1;
end
[U,S,V] = mySVD(X);
% S = S.^2;
% s = sort(diag(S),'descend');
% sumvalue = sum(s) * pcaRatio ;
% sumNow = 0;
% for i = 1:size(s)
%     sumNow = sumNow + s(i);
%     if sumNow >= sumvalue
%         break;
%     end
% end
% if i <= final_dimension 
%     i = final_dimension*2;
% else
%     if i > final_dimension +20
%         i = final_dimension + 20;
%     end
% end
i = final_dimension*times;
pca_vector = U(:,1:i)';
data = pca_vector * X;
end