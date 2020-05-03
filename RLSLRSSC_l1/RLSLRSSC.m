function [C, history] = RLSLRSSC(Y, label_truth, r, lambda4, lambda,lambda1, lambda2, lambda3, ITERMAX, error_type)
% it solves the following problem
% min lambda||PY-AC-E||^2_F + lambda1||C||_* + lambda2||C||_1 +
% lambda3||Y-P'PY||^2_F+lambda4||E||^2_F
% subject to: C^T*1_n=1_n
%-----------------------------------------
%P: d*D  Y:D*N   A:r*k   C:k*N
num_label = length(unique(label_truth));

ADAPTIVE_ON = 1;
[d n] = size(Y);

%initial the laplace matrix
option.K = 7; %4 - 7
option.k = 4; %can be tuned
L = initial_laplace(Y,option);

%use pca to reduce the dimension before initial P
[pca,Y] = myPCA(Y,r,1.2); %1.2-2

%intial the P
P0 = RLSLRSSC_initial_P(Y, L, r);

P = P0;
%D0 = dictnormalize(randn(r,k));
%D0 = DCT_D(Y,k);
%[D0,S] = Robust_pca(Y,25);
%D0 = dictnormalize(D0);
%D0 = Initializa_dictionary(Y,k,20);

%initial the noise matrix
E = zeros(size(P*Y));

%%%%% save the history acc
ACC = zeros(ITERMAX,1);
NMI = zeros(ITERMAX,1);
% lambda4 = 0.1;
for i = 1 : ITERMAX
    [C,History] = ADMM_RLSLRSSC( P*Y, P*Y-E,lambda, lambda1, lambda2, ADAPTIVE_ON);%step1
    %[C,H] = my_ALM_noisyLRSSC( A, P*Y,lambda, lambda1, lambda2, ADAPTIVE_ON);%step1
    P_new = RLSLRSSC_dimreduce(P, Y, C, E, L, r, lambda, lambda3);%step2
    
    c = 1;
    P = (1-c) * P + c * P_new;
    

%     A =  compute_dictionary_ODL(P*Y-E, A, C );%step3
    %PA = (P*Y)*pinv(C);
%     PA = my_ksvd_once( P*Y, P*A, C );
%     PA = compute_dictionary_MOD( P*Y, C );
    %A = metaface(P,Y,A,C);
%     E = lambda*(P*Y-P*Y*C)/(lambda+lambda4);
    E = update_E(P,Y,C,lambda,lambda4, error_type);
    
%     A = UpdateD_ALM(P*Y,A,C,lambda,lambda4);
    
    [label_cluster,NM,AC] = get_NMI_AC(label_truth, num_label, C);
    disp(['Iter : ' num2str(i) ', the NMI is : ' num2str(NM) ', the AC is : ' num2str(AC)]);
    cost = lambda*norm(P*Y-P*Y*C-E,'fro')^2 + lambda1*nuclear_norm(C)+ lambda2*norm(C,1) + lambda3*trace(P*Y*L*Y'*P')+lambda4 * E_norm(E,error_type);
    disp(['Iter : ' num2str(i) ', the cost is : ' num2str(cost)]);
    
    %save the history acc
    ACC(i) = AC;
    NMI(i) = NM;
end
history.ACC = ACC;
history.NMI = NMI;
return;



function [label_cluster, NMI, AC] = get_NMI_AC(label_truth, num_label, C)
W = BuildAdjacency(C);
% W = abs(C) + abs(C');
label_cluster = SpectralClustering(W, num_label);
[NMI, AC] = compute_NMI(label_truth, label_cluster);
return;