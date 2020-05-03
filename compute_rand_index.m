function rand_index=compute_rand_index(truth_label,cluster_label)
n=length(truth_label);
correct=0;
for i=1:n
    for j=(i+1):n
        if (truth_label(i)==truth_label(j))
            if(cluster_label(i)==cluster_label(j))
                correct=correct+1;
            end
        else
            if(cluster_label(i)~=cluster_label(j))
                correct=correct+1;
            end
        end
    end
end
rand_index=(2*correct)/(n*(n-1));
return;