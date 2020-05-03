function [RandIndex MIhat AC]=compute_RandIndex_NMI_and_AC(gnd,res)
RandIndex = compute_rand_index(gnd, res);
[N D] = size(gnd);
if (N == 1)
    gnd = gnd';
end
[N D] = size(res);
if (N == 1)
    res = res';
end
res = bestMap(gnd,res);
%=============  evaluate AC: accuracy ==============
AC = length(find(gnd == res))/length(gnd);
%=============  evaluate MIhat: nomalized mutual information =================
MIhat = MutualInfo(gnd,res);
return;

% function MIhat=compute_NMI(fea, gnd)
% fea = rand(50,70);
% gnd = [ones(10,1);ones(15,1)*2;ones(10,1)*3;ones(15,1)*4];
% 
% res = kmeans(fea,4);
% 
% res = bestMap(gnd,res);
% %=============  evaluate AC: accuracy ==============
% AC = length(find(gnd == res))/length(gnd);
% %=============  evaluate MIhat: nomalized mutual information =================
% MIhat = MutualInfo(gnd,res);
% return;