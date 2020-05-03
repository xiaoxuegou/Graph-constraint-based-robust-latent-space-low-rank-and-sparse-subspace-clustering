function [final_result] = litcount(filename)
fid = fopen(filename, 'rt');
RandIndex = [];
NMI = [];
AC = [];
while feof(fid) == 0
     line = fgetl(fid);
     effect_vec = strtok(line,'truth:');
     truth = str2num(effect_vec);
     line = fgetl(fid);
     effect_vec = strtok(line,'cluster:');
     cluster = str2num(effect_vec);
     
     
     [temp_RandIndex temp_NMI temp_AC] = compute_RandIndex_NMI_and_AC(truth, cluster);
     RandIndex = [RandIndex temp_RandIndex];
     NMI = [NMI temp_NMI];
     AC = [AC temp_AC];
end
fclose(fid);
final_result.RandIndex = batch_process(RandIndex);
final_result.NMI = batch_process(NMI);
final_result.AC = batch_process(AC);
return;

function [summary] = batch_process(result_array)
len = 5;
len_total = length(result_array);
num_fold = len_total/len;
col = zeros(len+2,1);
summary = [];
for i = 1:num_fold
    base = (i-1)*len;
    temp = result_array(base+1:base+5);
    avg_ = mean(temp);
    std_ = std(temp);
    col = [temp avg_ std_]';
    summary = [summary, col];
end
return;