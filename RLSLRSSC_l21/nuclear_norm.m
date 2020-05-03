function [nuclear] = nuclear_norm(M)
    [s,v,d] = svd(M);
    nuclear = sum(diag(v));
end