function [idx, idxS, idxC] = tcluster_components(S_time, K, num_clusters)
[~, ~, n3] = size(S_time);
Features = zeros(n3, K); 
for k = 1:K
    Features(:, k) = squeeze(S_time(k,k,:));
end

dist_temp = pdist(Features');
dist = squareform(dist_temp);
S = exp(-dist.^2);
issymmetric(S);
rng('default') 
idxS = spectralcluster(S,num_clusters,'Distance','precomputed','LaplacianNormalization','symmetric');

% 1. Calculate pairwise distances (returns a vector, not a square matrix)
dist_vec = pdist(Features'); 
% 2. Build the Hierarchical Tree using Ward's minimum variance method
Z = linkage(dist_vec, 'ward');
% 3. Cut the tree into exactly 3 clusters
idxC = cluster(Z, 'maxclust', num_clusters);


% Using K-means per your original script
idxK = kmeans(Features', num_clusters, 'Distance', 'sqeuclidean', 'Replicates', 10);

idx = idxK;
end