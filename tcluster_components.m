function idx = tcluster_components(S_time, K, num_clusters)
[~, ~, n3] = size(S_time);
Features = zeros(n3, K); 
for k = 1:K
    Features(:, k) = squeeze(S_time(k,k,:));
end
% Using K-means per your original script
idx = kmeans(Features', num_clusters, 'Distance', 'sqeuclidean', 'Replicates', 10);
end