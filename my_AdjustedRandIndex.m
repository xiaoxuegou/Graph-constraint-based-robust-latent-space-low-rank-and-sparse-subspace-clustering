function [ARI, RI] = my_AdjustedRandIndex(L1,L2)
% 编写代码的来源及出处：
% Vinh N X, Epps J, Bailey J. 
% Information Theoretic Measures for Clusterings Comparison: Variants, Properties, Normalization and Correction for Chance[M]. JMLR.org, 2010.

L1 = L1(:);
L2 = L2(:);

if size(L1) ~= size(L2)
    error('size(L1) must == size(L2)');
end

N = length(L1);
N00 = 0;%在L1中不属于同一类，在L2中不属于同一类
N01 = 0;%在L1中不属于同一类，在L2中属于同一类
N10 = 0;%在L1中属于同一类，在L2中不属于同一类
N11 = 0;%在L1中属于同一类，在L2中属于同一类

pairs_all = nchoosek(1:N,2);
a = pairs_all(:,1);
b = pairs_all(:,2);
len = length(a);

for t = 1:len
    index_i = a(t);
    index_j = b(t);
    if (L1(index_i)==L1(index_j)) %(1,?)
        if(L2(index_i)==L2(index_j))% (1,1)
            N11 = N11 + 1;
        else % (1,0)
            N10 = N10 + 1;
        end
    else % (0,?)
        if(L2(index_i)==L2(index_j))%  (0,1)
            N01 = N01 + 1;
        else% (0,0)
            N00 = N00 + 1;
        end
    end
end

RI = 2*(N00 + N11)/(N*(N-1));

ARI = 2*(N00*N11 - N01*N10)/((N00+N01)*(N01+N11) + (N00+N10)*(N10+N11));

return;