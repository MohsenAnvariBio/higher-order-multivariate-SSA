clear; close all; clc;

%% 0. Parameters & Data Loading (Prepared for Real ECG Data)
W = 1250;         % Window length (Set to 400 as per the paper's ECG experiment)
delta = 1;       % Step size 
num_clusters = 3; % 3 groups: Maternal ECG, Fetal ECG, and Noise

% =========================================================================
% REPLACE THIS BLOCK WITH YOUR REAL DATA LOADING
% Assuming you have a CSV or MAT file named 'real_ecg.csv' of size 2500 x 9
raw_data = load("FOETAL_ECG.mat"); 
raw_data = raw_data.FOETAL_ECG;
t = raw_data(:, 1)';         % Column 1 is time
data = raw_data(:, 2:9)';    % Columns 2-9 are the 8 channels. Transpose to (M x N)
% =========================================================================

% --- Temporary Mock Real Data (Remove when you load your real data) ---
N = 2500;  
t = 1:N;
M = 8;     % 8 Channels
% frequencies = linspace(1, 10, M)'; 
% data = sin(frequencies * t) + 0.2 * randn(M, N); % Mock signal
% ----------------------------------------------------------------------

%% 1. Step 1: Construct Trajectory Tensor
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

%% 2. Step 2: Tensor SVD (Algorithm 1 - Fourier Domain Decomposition)
[n1, n2, n3] = size(X);


% K = min(n1, n2); % 'Economy' rank limit
K = tubalrank(X); 
fprintf('Estimated Tubal Rank (K): %d\n', K);

X_fft = fft(X, [], 3); 

% Pre-allocate based on the economy size K
U_fft = zeros(n1, K, n3);   
S_fft = zeros(K, K, n3);    
V_fft = zeros(n2, K, n3);   

half_n3 = ceil((n3 + 1) / 2);
for i = 1:half_n3
    [U_i, S_i, V_i] = svd(X_fft(:,:,i), 'econ');
    U_fft(:,:,i) = U_i;
    S_fft(:,:,i) = S_i;
    V_fft(:,:,i) = V_i;
end

% Conjugate Symmetry for the second half
for j = (half_n3 + 1):n3
    U_fft(:,:,j) = conj(U_fft(:,:, n3 - j + 2));
    S_fft(:,:,j) = S_fft(:,:, n3 - j + 2);
    V_fft(:,:,j) = conj(V_fft(:,:, n3 - j + 2));
end

% Inverse FFT to get time-domain tensors
U = ifft(U_fft, [], 3);
S = ifft(S_fft, [], 3);
V = ifft(V_fft, [], 3);

%% 3. Step 3: Grouping via Spectral Clustering
% We must extract the "tubes" from the S tensor in the TIME domain (real part)
S_time = real(S);
Features = zeros(n3, K); 
for k = 1:K
    Features(:, k) = squeeze(S_time(k,k,:));
end

% Implement the exact Gaussian similarity matrix from the paper
dist_mat = squareform(pdist(Features')); % Distance between the K components
sigma = 1; 
A = exp(-dist_mat.^2 / (2 * sigma^2));

% Cluster components into 3 groups (Maternal, Fetal, Noise)
idx = spectralcluster(A, num_clusters, 'Distance', 'precomputed', ...
    'LaplacianNormalization', 'symmetric');
%% 4. Step 4: Reconstruction for ALL Clusters
num_clusters = 3;

% Pre-allocate a 3D matrix to hold [Cluster x Channel x Time]
reconstructed_data_all = zeros(num_clusters, M, N);

for c = 1:num_clusters
    % Create a mask for the current cluster in the loop
    component_mask = (idx == c);

    % 4a. Reconstruct the Tensor in the Fourier Domain
    X_rec_fft = zeros(size(X_fft));
    for i = 1:n3
        S_slice = S_fft(:,:,i);

        % Zero out the singular values NOT in our current cluster
        S_slice(~component_mask, ~component_mask) = 0;

        % Reconstruct the slice: X = U * S * V'
        X_rec_fft(:,:,i) = U_fft(:,:,i) * S_slice * V_fft(:,:,i)';
    end

    % 4b. Inverse FFT back to Time Domain
    X_rec = real(ifft(X_rec_fft, [], 3));

    % 4c. Block Diagonal Averaging for all channels in this cluster
    for m = 1:M
        reconstructed_data_all(c, m, :) = block_diagonal_averaging(X_rec(:,:,m), N, W, delta);
    end
end

%% 5. Visualization: 8 Figures (One per Channel)
% Display the first 1000 samples for better visibility of waveforms
view_range = 1:min(1000, N); 
colors = ['r', 'b', 'm']; % Colors to distinguish the 3 clusters

for m = 1:M
    % Create a new figure for each channel
    figure('Name', ['HO-MSSA: Channel ', num2str(m)], 'Color', 'w');

    % --- Subplot 1: Original Mixture ---
    subplot(num_clusters + 1, 1, 1);
    plot(t(view_range), data(m, view_range), 'k', 'LineWidth', 1); 
    title(['Channel ', num2str(m), ' - Original Mixture']);
    ylabel('Amplitude');
    grid on;
    set(gca, 'XTickLabel', []); % Hide X-ticks for top plots

    % --- Subplots 2 to 4: The 3 Extracted Clusters ---
    for c = 1:num_clusters
        subplot(num_clusters + 1, 1, c + 1);

        % Extract the 1D signal for this specific cluster and channel
        extracted_signal = squeeze(reconstructed_data_all(c, m, view_range));

        plot(t(view_range), extracted_signal, 'Color', colors(c), 'LineWidth', 1.2);
        title(['Cluster ', num2str(c), ' (Extracted Component)']);
        ylabel('Amplitude');
        grid on;

        % Only add the X-label to the very bottom plot
        if c == num_clusters
            xlabel('Time (Samples)');
        else
            set(gca, 'XTickLabel', []);
        end
    end

    % Adjust the window size so the subplots aren't cramped
    % The position shifts slightly for each m so the windows cascade cleanly
    set(gcf, 'Position', [50 + m*20, 50 + m*20, 700, 800]);
end

%% Function Definition: Block Diagonal Averaging
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