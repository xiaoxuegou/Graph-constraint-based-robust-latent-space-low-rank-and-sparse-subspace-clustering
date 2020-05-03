function NMI = my_MutualInfo(L1,L2)
%Vinh N X, Epps J, Bailey J. 
%Information Theoretic Measures for Clusterings Comparison: Variants, Properties, Normalization and Correction for Chance[M]. JMLR.org, 2010.

L1 = L1(:);
L2 = L2(:);

if size(L1) ~= size(L2)
    error('size(L1) must == size(L2)');
end

N = length(L1);

U = unique(L1);
num_U = length(U);
V = unique(L2);
num_V = length(V);

G = zeros(num_U, num_V);
for i=1:num_U
    for j=1:num_V
        G(i,j) = sum(L1 == U(i) & L2 == V(j));
    end
end

a = sum(G,2);%列向量，对每一行求和
b = sum(G,1); %行向量，对每一列求和

p_G = G/N;
p_a = a/N;
p_b = b/N;

H_U = sum(-p_a.*log(p_a)); %H(U)
H_V = sum(-p_b.*log(p_b)); % H(V)

% H_U_and_V = sum(sum(-p_G.*log(p_G))); % H(U,V)，算NMI_joint才用着
% ptemp = G./repmat(b, num_U, 1);
% H_U_frac_V = sum(-p_G(:).*log(ptemp(:))); %H(U/V)，算NMI用不着

ppp = N*G./(a*b);
ppp(abs(ppp) < 1e-12) = 1;%排除出现NaN的情况
I_U_and_V = sum(p_G(:).*log(ppp(:)));%I(U,V)


MI = I_U_and_V;

NMI_max = I_U_and_V/max(H_U, H_V);
% NMI_min = I_U_and_V/min(H_U, H_V);
% NMI_joint = I_U_and_V/H_U_and_V;
% NMI_sum = 2*I_U_and_V/(H_U + H_V);
% NMI_sqrt = I_U_and_V/sqrt(H_U*H_V);

NMI = NMI_max;

return;