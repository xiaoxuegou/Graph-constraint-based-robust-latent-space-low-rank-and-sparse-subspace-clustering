function [ output_args ] =  run_rlslrssc_l21_on_figure(db_name, process_mode,option)
dir_base = '~/RLSLRSSC/';
algorithm_name = 'RLSLRSSC_L21';
if nargin < 2
db_name = 'YaleB';
process_mode = 'Gaussian';
end
db = [db_name  '/'  process_mode];
dir_result_base = [dir_base '/experiment/result/figure/' db ];
if(~isdir(dir_result_base))
    mkdir(dir_result_base);
end
% dir_base = '~/pujulei';

k_this = [0:1:0];
%k_this = [1:1:10];
len_result = 1;
if nargin == 2
    for reducedDim=[100]
        for lambda = [0.01,0.001:0.001:0.005]
            for lambda1=[0.1]
                for lambda2= [lambda1/10]
                    for lambda3=[0.01]
                        for lambda4 = [0.1,0.05,0.02,0.01,0.005,0.002,0.001]
                            for ITERMAX=[5]
                                for error_type = [21]
                                    run_throughout(dir_base, dir_result_base, db, algorithm_name, len_result, k_this, reducedDim, lambda, lambda1, lambda2, lambda3, ITERMAX, lambda4, error_type);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
if nargin == 3
    error_type = 21;
    k_this = option.k;
    len_result = option.len_result;
    reducedDim = option.reducedDim;
    lambda = option.lambda;
    lambda1 = option.lambda1;
    lambda2 = option.lambda2;
    lambda3 = option.lambda3;
    lambda4 = option.lambda4;
    ITERMAX = option.ITERMAX;
    run_throughout(dir_base, dir_result_base, db, algorithm_name, len_result, k_this, reducedDim, lambda, lambda1, lambda2, lambda3, ITERMAX, lambda4, error_type);
end
return;


function[] = run_throughout(dir_base, dir_result_base, db, algorithm_name, len_result, k_this, reducedDim, lambda,lambda1, lambda2, lambda3, ITERMAX, lambda4, error_type)
for i = k_this    
dir_result = [dir_result_base '/K=' num2str(i) '/' algorithm_name '/r=' num2str(reducedDim) ',lambda= ' num2str(lambda)  ',lambda1= ' num2str(lambda1) ',lambda2= ' num2str(lambda2) ',lambda3=' num2str(lambda3) ',lambda4=' num2str(lambda4) ',ITERMAX=' num2str(ITERMAX)];
if(~isdir(dir_result))
    mkdir(dir_result);
end

result_record_name = [dir_result '/meta_result.txt'];
result_record_fid = fopen(result_record_name, 'w+');
prt2file(result_record_fid, ['num_trial:' num2str(len_result)]);


dir_data = [dir_base '/experiment/data/figure/' db '/K=' num2str(i)];
run_once(result_record_fid ,reducedDim, lambda,lambda1, lambda2, lambda3, ITERMAX, lambda4,error_type, dir_data, dir_result, len_result); 
end

return;

function [] = run_once(result_record_fid ,reducedDim, lambda, lambda1, lambda2, lambda3, ITERMAX, lambda4, error_type, dir_data, dir_result, len_result)
flag_run = true;
% max_iter = 400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% record_name = make_record_name(db_name, dir_result, reducedDim, lambda1, lambda2, lambda3, ITERMAX, atomSizePercent); 
record_name = [dir_result '/' 'result.txt'];
fid = fopen(record_name, 'w+');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AC = zeros(len_result,1);%
NMI = zeros(len_result,1);%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idx =1:1:len_result
    clear data_train data_test label_train label_test;
%     fileName=[ 'cross_validation_' num2str(current_fold)];
    fileName=[ 'cross_validation_' num2str(1)];
    load([dir_data '/' fileName]);
    Y=data_train;
    label_truth=label_train;
    num_label=length(unique(label_train));
    label_cluster=zeros(1, length(label_train));
    
    if (flag_run == true)
        [D N] = size(data_train);
        [C history] = RLSLRSSC(Y, label_train, reducedDim, lambda4, lambda, lambda1, lambda2, lambda3, ITERMAX, error_type);
        [label_cluster,NMI(idx), AC(idx)] = getSegmentation(label_truth, num_label, C);
    end
    
    prt2file(result_record_fid, ['truth:' num2str(convert_to_row(label_truth))]);
    prt2file(result_record_fid, ['cluster:' num2str(convert_to_row(label_cluster))]);
    
    str_time=get_str_time();
    
    prt(fid, ['AC:' num2str(AC(idx)) ', NMI:' num2str(NMI(idx))]);
    prt(fid, ['All iterative ACC:' num2str(convert_to_row(history.ACC))]);
    prt(fid, ['All iterative NMI:' num2str(convert_to_row(history.NMI))]);
    prt(fid, ['************************************************************************************']);
end
prt(fid, ['****************************' ' summary ' '****************************']);
prt(fid, ['average_AC:' num2str(mean(AC)) ',std_AC:' num2str(std(AC))]);
prt(fid, ['average_NMI:' num2str(mean(NMI)) ',std_NMI:' num2str(std(NMI))]);
fclose(fid);
return;

function [label_cluster, NMI, ACC] = getSegmentation(label_truth, num_label, C)
W = BuildAdjacency(C);
label_cluster = SpectralClustering(W, num_label);
[NMI ,ACC] = compute_NMI(label_truth, label_cluster);
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
