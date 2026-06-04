clear; close all; clc;


%% 0. Parameters & Data Loading
% sin
% % --- Temporary Mock Real Data (Remove when you load your real data) ---
% N = 2500;  
% % t = 1:N;
% M = 5;     % 5 Channels
% % frequencies = linspace(1, 10, M)'; 
% % data = sin(frequencies * t) + 0.2 * randn(M, N); % Mock signal
% % ----------------------------------------------------------------------

% ECG
% W = 400;         
% delta = 1;       
% num_clusters = 3; 
% N = 800; 
% raw_data = load("FOETAL_ECG.mat"); 
% raw_data = raw_data.FOETAL_ECG;
% % raw_data(1:N, 3) = -raw_data(1:N, 3);
% t = raw_data(1:N, 1)';         
% data = raw_data(1:N, 2:6)';    
% M = 5; 
% % CRITICAL FIX 1: Normalize globally, NOT channel-by-channel.
% % This preserves the physical spatial mixing matrix of the electrodes.
% data = data ./ max(abs(data(:)));

W = 1500;         
delta = 1;       
num_clusters = 2; 
fs = 200; 
load("Contaminated_Data.mat"); 
data = sim10_con([1, 2, 3, 4, 11, 12], 1:3000); %EEG channels (Fp1, Fp2, F3, F4, F7, F8)
[M, N] = size(data); 
t = 1: 1/fs: N;         
% CRITICAL FIX 1: Normalize globally, NOT channel-by-channel.
% This preserves the physical spatial mixing matrix of the electrodes.
data = data ./ max(abs(data(:)));

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
% K = 92; 
X_fft = fft(X, [], 3); 

% Pre-allocate based on the economy size K
U_fft = zeros(n1, K, n3);   
S_fft = zeros(K, K, n3);    
V_fft = zeros(n2, K, n3);   

half_n3 = ceil((n3 + 1) / 2);
for i = 1:half_n3
    [U_i, S_i, V_i] = svd(X_fft(:,:,i), 'econ');
% CRITICAL: Explicitly slice the outputs to keep only K components
    U_fft(:,:,i) = U_i(:, 1:K);
    S_fft(:,:,i) = S_i(1:K, 1:K);
    V_fft(:,:,i) = V_i(:, 1:K);
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

% [U_tsvd,S_tsvd,V_tsvd] = tsvd(X,'full');

%% 3. Step 3: Grouping via Spectral Clustering
S_time = real(S);
Features = zeros(n3, K); 
for k = 1:K
    Features(:, k) = squeeze(S_time(k,k,:));
end

dist_temp = pdist(Features');
dist = squareform(dist_temp);
S = exp(-dist.^2);
issymmetric(S);
rng('default') 
idxS1 = spectralcluster(S,num_clusters,'Distance','precomputed','LaplacianNormalization','symmetric');

% Features = Features ./ max(abs(Features(:)));

idxS = spectralcluster(Features, num_clusters);

idxK = kmeans(Features', num_clusters, 'Distance', 'sqeuclidean', 'Replicates', 10);

% 1. Calculate pairwise distances (returns a vector, not a square matrix)
dist_vec = pdist(Features'); 
% 2. Build the Hierarchical Tree using Ward's minimum variance method
Z = linkage(dist_vec, 'ward');
% 3. Cut the tree into exactly 3 clusters
idxC = cluster(Z, 'maxclust', num_clusters);

idx = idxK;
%% 4. Step 4: Reconstruction for ALL Clusters
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
% Display the first 800 samples for better visibility of waveforms
view_range = 1:min(4000, N); 
colors = ['r', 'b', 'm']; % Colors to distinguish the 3 clusters

for m = 1:M
    % Create a new figure for each channel
    figure('Name', ['HO-MSSA: Channel ', num2str(m)], ...
             'Color', 'w');
    % movegui(fig, 'center');

    % --- Subplot 1: Original Mixture ---
    subplot(num_clusters + 1, 1, 1);
    plot(t(view_range), data(m, view_range), 'k', 'LineWidth', 1); 
    title(['Channel ', num2str(m), ' - Original Mixture']);
    ylabel('Amplitude');
    y_limits = ylim;
    % ylim([-1, 1])
    grid on;
    set(gca, 'XTickLabel', []); % Hide X-ticks for top plots

    % --- Subplots 2 to 4: The 3 Extracted Clusters ---
    for c = 1:num_clusters
        subplot(num_clusters + 1, 1, c + 1);

        % Extract the 1D signal for this specific cluster and channel
        extracted_signal = squeeze(reconstructed_data_all(c, m, view_range));

        plot(t(view_range), extracted_signal, 'Color', colors(c), 'LineWidth', 1.2);
        title(['Cluster ', num2str(c)]);
        ylabel('Amplitude');
        ylim(y_limits);
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
    set(gcf, 'Position', [50 + m*20, 50 + m*20, 600, 300]);
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