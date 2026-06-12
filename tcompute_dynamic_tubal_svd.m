function [U_fft, S_fft, V_fft, S_time, K, s_k_norm, VR_k] = compute_dynamic_tubal_svd(X, threshold)
[n1, n2, n3] = size(X);
I_dim = min(n1, n2); 
half_n3 = ceil((n3 + 1) / 2);
X_fft = fft(X, [], 3); 

% Compute norms to find K
S_full_fft = zeros(I_dim, I_dim, n3);
for i = 1:half_n3
    [~, S_i, ~] = svd(X_fft(:,:,i), 'econ');
    S_full_fft(1:size(S_i,1), 1:size(S_i,2), i) = S_i; 
end
for j = (half_n3 + 1):n3
    S_full_fft(:,:,j) = S_full_fft(:,:, n3 - j + 2);
end

S_full_time = real(ifft(S_full_fft, [], 3));
s_k_norm = zeros(I_dim, 1);
for k = 1:I_dim
    s_k_norm(k) = norm(squeeze(S_full_time(k, k, :)), 2); 
end

s_k_norm_sq = s_k_norm.^2;
VR_k = cumsum(s_k_norm_sq) / sum(s_k_norm_sq);
K = find(VR_k >= threshold, 1);
if isempty(K), K = I_dim; end

% Truncated SVD using calculated K
U_fft = zeros(n1, K, n3);   
S_fft = zeros(K, K, n3);    
V_fft = zeros(n2, K, n3);   
for i = 1:half_n3
    [U_i, S_i, V_i] = svd(X_fft(:,:,i), 'econ');
    U_fft(:,:,i) = U_i(:, 1:K);
    S_fft(:,:,i) = S_i(1:K, 1:K);
    V_fft(:,:,i) = V_i(:, 1:K);
end
for j = (half_n3 + 1):n3
    U_fft(:,:,j) = conj(U_fft(:,:, n3 - j + 2));
    S_fft(:,:,j) = S_fft(:,:, n3 - j + 2);
    V_fft(:,:,j) = conj(V_fft(:,:, n3 - j + 2));
end

% Time domain S for clustering
S_time = real(ifft(S_fft, [], 3));
end
