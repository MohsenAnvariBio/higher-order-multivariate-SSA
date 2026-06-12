function [s_k_norm, VR_k, K] = tcompute_tubal_metrics(X, threshold)
% Computes the singular tube L2 norms, Variance Ratio, and dynamic K

[n1, n2, n3] = size(X);
I_dim = min(n1, n2); 
half_n3 = ceil((n3 + 1) / 2);

% FFT along 3rd dimension
X_fft = fft(X, [], 3); 
S_full_fft = zeros(I_dim, I_dim, n3);

% Calculate full economy S
for i = 1:half_n3
    [~, S_i, ~] = svd(X_fft(:,:,i), 'econ');
    S_full_fft(1:size(S_i,1), 1:size(S_i,2), i) = S_i; 
end

% Conjugate symmetry for the second half
for j = (half_n3 + 1):n3
    S_full_fft(:,:,j) = S_full_fft(:,:, n3 - j + 2);
end

% Bring S back to time domain
S_full_time = real(ifft(S_full_fft, [], 3));

% Pre-allocate and calculate L2 norm of each singular tube ||s_k||_2
s_k_norm = zeros(I_dim, 1);
for k = 1:I_dim
    tube_k = squeeze(S_full_time(k, k, :));
    s_k_norm(k) = norm(tube_k, 2); 
end

% Calculate Variance Ratio VR(k)
s_k_norm_sq = s_k_norm.^2;
total_variance = sum(s_k_norm_sq);
VR_k = cumsum(s_k_norm_sq) / total_variance;

% Determine K based on threshold
K = find(VR_k >= threshold, 1);
if isempty(K)
    K = I_dim; 
end
end