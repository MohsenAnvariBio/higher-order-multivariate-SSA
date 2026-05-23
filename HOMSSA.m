clear; close all; clc;

%% 0. Parameters & Data Generation (Modified for visibility)
W = 200;         % Window length
delta = 1;      % Step size
M = 5;          % Number of channels (e.g., 5 EEG sensors)
N = 5000;        % Length of signal
t = 1:N;

% Generate a structured signal: Sine wave (Common to all) + Random Noise
% pure_signal = sin(0.1* t);
% data = repmat(pure_signal, M, 1) + 0.15 * randn(M, N) + 0.25 * randn(M, N); 
frequencies = [10; 12; 15; 17; 20]; % Column vector (5x1)
pure_signals = sin(frequencies * t);
data = pure_signals + 0.25 * randn(M, N);

%% 1. Step 2: Construct Trajectory Tensor
% Dimensions: [WindowLength x NumColumns x NumChannels]
num_cols = floor((N - W) / delta) + 1;
X = zeros(W, num_cols, M);

for m = 1:M
    % Using the paper's TDE logic
    for j = 1:num_cols
        start_idx = (j-1)*delta + 1;
        X(:, j, m) = data(m, start_idx : start_idx + W - 1);
    end
end

%% 2. Algorithm 1: t-SVD (Fourier Domain Decomposition)
[n1, n2, n3] = size(X);
K = min(n1, n2); % This is the 'Economy' rank (should be 20 in my case)

X_fft = fft(X, [], 3); 

% Pre-allocate based on the economy size K
U_fft = zeros(n1, K, n3);   
S_fft = zeros(K, K, n3);    
V_fft = zeros(n2, K, n3);   

% Step 2: SVD on the first half of slices
half_n3 = ceil((n3 + 1) / 2);
for i = 1:half_n3
    % 'econ' returns U(n1 x K), S(K x K), V(n2 x K)
    [U_i, S_i, V_i] = svd(X_fft(:,:,i), 'econ');
    
    U_fft(:,:,i) = U_i;
    S_fft(:,:,i) = S_i;
    V_fft(:,:,i) = V_i;
end

% Step 3: Fill the second half (Conjugate Symmetry)
for j = (half_n3 + 1):n3
    U_fft(:,:,j) = conj(U_fft(:,:, n3 - j + 2));
    S_fft(:,:,j) = S_fft(:,:, n3 - j + 2);
    V_fft(:,:,j) = conj(V_fft(:,:, n3 - j + 2));
end

% Step 4: Inverse FFT
U = ifft(U_fft, [], 3);
S = ifft(S_fft, [], 3);
V = ifft(V_fft, [], 3);


%% 3. Step 3: Grouping via Spectral Clustering
% We extract the "tubes" from the S tensor in the frequency domain
K = min(n1, n2); 
Features = zeros(n3, K); 
for k = 1:K
    Features(:, k) = squeeze(S_fft(k,k,:));
end

% Clustering components into 2 groups (Signal vs Noise)
% idx tells us which component (1 to K) belongs to which group
% dist = squareform(pdist(S));
% A = exp(-dist.^2);
% idx = spectralcluster(A,2,'Distance','precomputed','LaplacianNormalization','symmetric');
idx = spectralcluster(real(Features'), 2); 
gscatter(Features(1,:),Features(5,:),idx);
% Identify which cluster is the "Signal" (usually the one with more energy)
% We'll assume Cluster 1 is signal for this example
target_cluster = 2; 
component_mask = (idx == target_cluster);

%% 4. Step 4: Reconstruction
% 4a. Reconstruct the Tensor in the Fourier Domain
% We only keep the SVD components that belong to our target cluster
X_rec_fft = zeros(size(X_fft));
for i = 1:n3
    % Pull out the S matrix for this slice
    S_slice = S_fft(:,:,i);
    % Zero out the singular values NOT in our cluster
    % component_mask is 1xK, we apply it to the diagonal
    S_slice(~component_mask, ~component_mask) = 0;
    
    % Reconstruct the slice: X = U * S * V'
    X_rec_fft(:,:,i) = U_fft(:,:,i) * S_slice * V_fft(:,:,i)';
end

% 4b. Inverse FFT back to Time Domain
X_rec = real(ifft(X_rec_fft, [], 3));

% 4c. Block Diagonal Averaging (Final Signal Recovery)
reconstructed_data = zeros(M, N);
for m = 1:M
    reconstructed_data(m, :) = block_diagonal_averaging(X_rec(:,:,m), N, W, delta);
end

% %% 5. Visualization
% figure('Color', 'w');
% subplot(2,1,1);
% plot(t(1:1000), data(1,1:1000), 'Color', [0.7 0.1 0.7]); hold on;
% plot(t(1:1000), pure_signal(1,1:1000), 'b--', 'LineWidth', 0.1);
% title('Original Noisy Channel 1 vs Ground Truth');
% xlabel('Time Index (Samples)');
% ylabel('Amplitude');
% legend('Noisy Signal', 'Pure Sine');
% 
% subplot(2,1,2);
% plot(t(1:1000), reconstructed_data(1,1:1000), 'b', 'LineWidth', 1.5);
% title('HO-MSSA Reconstructed Channel 1');
% xlabel('Time Index (Samples)');
% ylabel('Amplitude');
% legend('Cleaned Signal');
%% 5. Visualization: All 5 Channels
figure('Name', 'HO-MSSA: Multi-Channel Reconstruction', 'Color', 'w');

% Define a small segment to view clearly
view_range = 1:200; 

for m = 1:M
    % Create a subplot for each channel
    subplot(M, 1, m);

    % Plot Noisy vs Reconstructed
    plot(t(view_range), data(m, view_range), 'Color', [0.8 0.8 0.8], 'DisplayName', 'Noisy'); 
    hold on;
    plot(t(view_range), 1*reconstructed_data(m, view_range), 'r', 'LineWidth', 0.5, 'DisplayName', 'Cleaned');
    plot(t(view_range), 1*pure_signals(m, view_range), 'g', 'LineWidth', 0.2, 'DisplayName', 'Cleaned');

    % Formatting each subplot
    ylabel(['Ch ', num2str(m)]);
    grid on;

    % Only add Title to the top plot
    if m == 1
        title('HO-MSSA Results: All 5 Channels (Noisy vs. Reconstructed)');
        legend('Location', 'northeast');
    end

    % Only add X-label to the bottom plot to save space
    if m == M
        xlabel('Time Index (Samples)');
    else
        set(gca, 'XTickLabel', []); % Remove X-ticks for middle plots
    end
end

% Optional: Adjust figure size for better vertical visibility
set(gcf, 'Position', [100, 100, 800, 900]);

%% Function Definition
function x_final = block_diagonal_averaging(X_slice, N, W, delta)
    [~, NumCols] = size(X_slice);
    x_sum = zeros(1, N);
    counts = zeros(1, N);
    for j = 1:NumCols
        start_idx = (j-1) * delta + 1;
        end_idx = start_idx + W - 1;
        if end_idx > N, end_idx = N; end
        actual_len = end_idx - start_idx + 1;
        x_sum(start_idx:end_idx) = x_sum(start_idx:end_idx) + X_slice(1:actual_len, j)';
        counts(start_idx:end_idx) = counts(start_idx:end_idx) + 1;
    end
    x_final = x_sum ./ counts;
end