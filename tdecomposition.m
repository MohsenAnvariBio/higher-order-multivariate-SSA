function [U_fft, S_fft, V_fft, S_time, K, s_k_norm, VR_k] = tdecomposition(X, variance_threshold, method)
% TCOMPUTE_DYNAMIC_TUBAL_SVD Computes t-SVD, Tucker, or SSA decompositions
% Inputs:
%   X         : 3D tensor of size n1 x n2 x n3
%   threshold : Variance threshold to determine truncation rank K
%   method    : (Optional) 'tsvd' (default), 'tucker', or 'ssa'

    % Default to original t-SVD if no method is specified
    if nargin < 3
        method = 'tsvd';
    end

    [n1, n2, n3] = size(X);
    I_dim = min(n1, n2); 
    half_n3 = ceil((n3 + 1) / 2);
    X_fft = fft(X, [], 3); 

    % =========================================================================
    % 1. Compute norms to find K (Kept globally identical for all methods)
    % =========================================================================
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
    K = find(VR_k >= variance_threshold, 1);
    if isempty(K), K = I_dim; end

    % Initialize outputs
    U_fft = zeros(n1, K, n3);   
    S_fft = zeros(K, K, n3);    
    V_fft = zeros(n2, K, n3);   

    % =========================================================================
    % 2. Execute chosen Decomposition Method
    % =========================================================================
    switch lower(method)
        
        case 'tsvd'
            % --- ORIGINAL T-SVD (Frequency Domain) ---
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
            S_time = real(ifft(S_fft, [], 3));

        case 'tucker'
            % --- TUCKER-BASED DECOMPOSITION (HOSVD mapping to tubal) ---
            % Unfold mode 1 to get spatial U basis
            X_unfold1 = reshape(X, n1, n2*n3);
            [U1, ~, ~] = svd(X_unfold1, 'econ');
            U1 = U1(:, 1:K);
            
            % Unfold mode 2 to get spatial V basis
            X_unfold2 = reshape(permute(X, [2, 1, 3]), n2, n1*n3);
            [U2, ~, ~] = svd(X_unfold2, 'econ');
            U2 = U2(:, 1:K);
            
            % Map to tubal time-domain format
            U_time = zeros(n1, K, n3);
            V_time = zeros(n2, K, n3);
            S_time = zeros(K, K, n3);
            
            % Concentrate spatial bases at t=1 for tubal convolution logic
            U_time(:,:,1) = U1; 
            V_time(:,:,1) = U2; 
            
            for i = 1:n3
                % The tubal Core tensor S
                S_time(:,:,i) = U1' * X(:,:,i) * U2;
            end
            
            % Transform to Fourier domain to match output signature
            U_fft = fft(U_time, [], 3);
            V_fft = fft(V_time, [], 3);
            S_fft = fft(S_time, [], 3);

        case 'ssa'
            % --- SSA-BASED DECOMPOSITION (Time-Domain Slice-wise SVD) ---
            U_time = zeros(n1, K, n3);
            V_time = zeros(n2, K, n3);
            S_time = zeros(K, K, n3);
            
            for i = 1:n3
                [U_i, S_i, V_i] = svd(X(:,:,i), 'econ');
                U_time(:,:,i) = U_i(:, 1:K);
                S_time(:,:,i) = S_i(1:K, 1:K);
                V_time(:,:,i) = V_i(:, 1:K);
            end
            
            % Transform to Fourier domain to match output signature
            U_fft = fft(U_time, [], 3);
            S_fft = fft(S_time, [], 3);
            V_fft = fft(V_time, [], 3);
            
        otherwise
            error('Invalid method. Please choose ''tsvd'', ''tucker'', or ''ssa''.');
    end
end