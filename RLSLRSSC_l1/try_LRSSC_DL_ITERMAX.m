function [ret] = try_LRSSC_DL_ITERMAX()
dir_base = 'F:/Æ×¾ÛÀà';
algorithm_name = 'my_algorithm';
db_name = 'COIL20';
data_type = 'figure';
num_fold = 1;
ITERMAX  = 10;
dir_result_base = [dir_base '/experiment/result/' data_type '/' db_name  '/' algorithm_name '/%time=' get_str_time()];
if(~isdir(dir_result_base))
    mkdir(dir_result_base);
end
result_record_name = [dir_result_base '/meta_result.txt'];
result_record_fid = fopen(result_record_name, 'w+');
prt2file(result_record_fid, ['num_trial:' num2str(ITERMAX+1)]);
record_name = [dir_result_base '/itermax.txt'];
fid = fopen(record_name, 'w+');
for cur_fold = 1:num_fold
    str_data = [dir_base  '/experiment/data/' data_type '/' db_name '/K=12/cross_validation_' num2str(cur_fold)];
%     str_data = [dir_base  '/experiment/data/' data_type '/' db_name '/K=8/cross_validation_1'];
    load(str_data);
    prt(fid, ['**********************validation_' num2str(cur_fold) '*****************************']);
    ret = try_once(result_record_fid, fid, data_train, label_train, ITERMAX);
end
return;

function [ret] = try_once(result_record_fid, fid, data_train, label_train, ITERMAX)
flag =  true;
Y = data_train;
label_truth = label_train;
num_label = length(unique(label_truth));
[d n] =size(Y);
r = 100;
lambda1 = 1; %low-rank
lambda2 = 0.01;%sparse
lambda3 = 1;
k = ceil(0.5*n);
label_cluster = zeros(1, n);
ADAPTIVE_ON = 1;
max_iter = 100;

if (flag == true)
    D0 = dictnormalize(rand(d,k));
    C0 = my_ALM_noisyLRSSC( D0, Y, lambda1, lambda2, ADAPTIVE_ON, max_iter);
    P0 = 0;
    
    A = D0;
    C = C0;
    P = P0;
    [tempAC, label_cluster] = wrap_func(fid, label_truth, num_label, C, 0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
prt2file(result_record_fid, ['truth:' num2str(convert_to_row(label_truth))]);
prt2file(result_record_fid, ['cluster:' num2str(convert_to_row(label_cluster))]);

for i = 1 : ITERMAX
    if (flag == true)
%         W = pinv(Y) * A *  C;
%         P = my_dimreduce(Y, W, r, 0.5, lambda3);%step1
%         PA =  compute_dictionary_ODL( P*Y, P*A, C );%step2
%         C = my_ALM_noisyLRSSC( PA, P*Y, lambda1, lambda2, ADAPTIVE_ON, max_iter);%step3
%         A = pinv(P)*PA;
        
        W = A *  C;
        A =  compute_dictionary_ODL( Y, A, C );%step2
        C = my_ALM_noisyLRSSC( A, Y, lambda1, lambda2, ADAPTIVE_ON, max_iter);%step3
        
        [tempAC, label_cluster] = wrap_func(fid, label_truth, num_label, C, i);
    end

    prt2file(result_record_fid, ['truth:' num2str(convert_to_row(label_truth))]);
    prt2file(result_record_fid, ['cluster:' num2str(convert_to_row(label_cluster))]);
end
ret = A*C;
return;

function [tempAC, label_cluster] = wrap_func(fid, label_truth, num_label, C, iter)
[tempAC, label_cluster] = special_func(fid, '|C^T*C|', label_truth, num_label, C, iter);
[tempAC, label_cluster] = special_func(fid, '|C^T|*|C|', label_truth, num_label, C, iter);
[tempAC, label_cluster] = special_func(fid, 'C^T*C', label_truth, num_label, C, iter);
return;

function [tempAC, label_cluster] = special_func(fid, name_method, label_truth, num_label, C, iter)
    if (strcmp(name_method, '|C^T*C|') == 1)
        W = abs(C'*C);
        [label_cluster tempAC] = getAC(label_truth, num_label, W);
        str = ['AC_' name_method ':' num2str(tempAC) '*******' 'ITERMAX=' num2str(iter)];
        prt(fid,str);
    elseif (strcmp(name_method, '|C^T|*|C|') == 1)
        W = abs(C')*abs(C);
        [label_cluster tempAC] = getAC(label_truth, num_label, W);
        str = ['AC_' name_method ':' num2str(tempAC) '*******' 'ITERMAX=' num2str(iter)];
        prt(fid,str);
    else
        name_method = 'C^T*C';
        W = C'*C;
        [label_cluster tempAC] = getAC(label_truth, num_label, W);
        str = ['AC_' name_method ':' num2str(tempAC) '*******' 'ITERMAX=' num2str(iter)];
        prt(fid,str);
    end
return;

function [label_cluster tempAC] = getAC(label_truth, num_label, W)
label_cluster = SpectralClustering(W, num_label);
[tempTrash tempAC] = compute_NMI(label_truth, label_cluster);
return;

function [] =  prt(fid, str)
disp(str);
fprintf(fid,'%s \r\n',str);
return;

function [] = prt2file(fid, str)
fprintf(fid,'%s \r\n',str);
return;

function [out] = convert_to_row(in)
[d n] = size(in);
if(d == 1)
    out = in;
else
    out = in';
end
return;

function [str_time] = get_str_time()
time_now=clock;
str_time = [num2str(time_now(1)) '_' num2str(time_now(2)) '_' num2str(time_now(3)) '_' num2str(time_now(4)) '_' num2str(time_now(5)) '_' num2str(time_now(6))];
return;