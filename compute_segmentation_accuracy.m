function segmentation_accuracy=compute_segmentation_accuracy(label_truth,label_cluster)
label_reshape = label_cluster;
label_cluster_unique = unique(label_cluster);
label_cluster_unique_reshape = [];
count = 0;
for i = 1:length(label_cluster_unique)
    t = label_cluster_unique(i);
    index = find(label_cluster==t);
    table = tabulate(label_truth(index));
    [maxCount,idx] = max(table(:,2));
    label_reshape(index) = 
    label_cluster_unique_reshape = [label_cluster_unique_reshape table(idx)];
end
segmentation_accuracy = length(find(label_truth == label_cluster_unique_reshape))/length(label_truth);
return;