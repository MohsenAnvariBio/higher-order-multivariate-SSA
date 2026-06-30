clear; close all; clc;

%% 0. User Input & Data Loading
% Prompt the user to select which dataset to run
patient_id = input('Enter the dataset number to process (1 to 27): ');

% Input validation
if isempty(patient_id) || patient_id < 1 || patient_id > 27
    fprintf('Invalid input. Defaulting to dataset 3.\n');
    patient_id = 3; 
end

% Parameters
W = 1500;         
delta = 1;       
num_clusters = 2; 
fs = 200; 
variance_threshold = 0.999; % 99.9% variance for dynamic K
channels_to_use = [1, 2, 3, 4, 11, 12]; % Fp1, Fp2, F3, F4, F7, F8
time_samples = 1:3000;

% Dynamically load the selected dataset
fprintf('Loading dataset sim%d_con...\n', patient_id);
all_data = load("Contaminated_Data.mat"); 
var_name = sprintf('sim%d_con', patient_id);

if ~isfield(all_data, var_name)
    error('Variable %s not found in Contaminated_Data.mat.', var_name);
end

raw_data = all_data.(var_name)(channels_to_use, time_samples);
[M, N] = size(raw_data); 
t = 1:1/fs:N; 

% CRITICAL FIX 1: Normalize globally, NOT channel-by-channel.
data = raw_data ;%./ max(abs(raw_data(:)));

%% 1. Step 1: Construct Trajectory Tensor
fprintf('Step 1: Constructing Trajectory Tensor...\n');
X = tbuild_trajectory_tensor(data, W, delta);

%% 2. Step 2: Tensor SVD & Figure 13
fprintf('Step 2: Performing Tubal SVD...\n');
[U_fft, S_fft, V_fft, S_time, K, s_k_norm, VR_k] = tcompute_dynamic_tubal_svd(X, variance_threshold);

fprintf('Calculated Tubal Rank K for patient %d (capturing %.2f%% variance): %d\n', patient_id, variance_threshold * 100, K);

% Plot Figure 13 for this specific dataset
tplot_figure_13_single(s_k_norm, VR_k, patient_id);

%% 3. Step 3: Grouping via Spectral Clustering / K-Means
fprintf('Step 3: Clustering Components...\n');
idx = tcluster_components(S_time, K, num_clusters);

%% 4. Step 4: Reconstruction for ALL Clusters
fprintf('Step 4: Reconstructing Signals...\n');
reconstructed_data_all = treconstruct_signals(U_fft, S_fft, V_fft, idx, num_clusters, M, N, W, delta);

%% 5. Visualization: Signal Waveforms
fprintf('Step 5: Generating Plots...\n');
view_range = 1:min(4000, N); 
colors = ['r', 'b', 'm', 'g']; % Expandable if num_clusters > 3

for m = 1:M
    figure('Name', ['HO-MSSA: Patient ', num2str(patient_id), ' - Channel ', num2str(m)], 'Color', 'w');
    
    % --- Subplot 1: Original Mixture ---
    subplot(num_clusters + 1, 1, 1);
    plot(t(view_range), data(m, view_range), 'k', 'LineWidth', 1); 
    title(['Channel ', num2str(m), ' - Original Mixture']);
    ylabel('Amplitude');
    y_limits = ylim;
    grid on;
    set(gca, 'XTickLabel', []); 
    
    % --- Subplots for Extracted Clusters ---
    for c = 1:num_clusters
        subplot(num_clusters + 1, 1, c + 1);
        extracted_signal = squeeze(reconstructed_data_all(c, m, view_range));
        plot(t(view_range), extracted_signal, 'Color', colors(mod(c-1, length(colors))+1), 'LineWidth', 1.2);
        title(['Cluster ', num2str(c)]);
        ylabel('Amplitude');
        ylim(y_limits);
        grid on;
        
        if c == num_clusters
            xlabel('Time (Samples)');
        else
            set(gca, 'XTickLabel', []);
        end
    end
    set(gcf, 'Position', [50 + m*20, 50 + m*20, 600, 300]);
end
fprintf('Processing Complete!\n');


%% HELPER FUNCTIONS
% function X = build_trajectory_tensor(data, W, delta)
%     [M, N] = size(data);
%     num_cols = floor((N - W) / delta) + 1;
%     X = zeros(W, num_cols, M);
%     for m = 1:M
%         for j = 1:num_cols
%             start_idx = (j-1)*delta + 1;
%             X(:, j, m) = data(m, start_idx : start_idx + W - 1);
%         end
%     end
% end
% 
% function [U_fft, S_fft, V_fft, S_time, K, s_k_norm, VR_k] = compute_dynamic_tubal_svd(X, threshold)
%     [n1, n2, n3] = size(X);
%     I_dim = min(n1, n2); 
%     half_n3 = ceil((n3 + 1) / 2);
%     X_fft = fft(X, [], 3); 
% 
%     % Compute norms to find K
%     S_full_fft = zeros(I_dim, I_dim, n3);
%     for i = 1:half_n3
%         [~, S_i, ~] = svd(X_fft(:,:,i), 'econ');
%         S_full_fft(1:size(S_i,1), 1:size(S_i,2), i) = S_i; 
%     end
%     for j = (half_n3 + 1):n3
%         S_full_fft(:,:,j) = S_full_fft(:,:, n3 - j + 2);
%     end
% 
%     S_full_time = real(ifft(S_full_fft, [], 3));
%     s_k_norm = zeros(I_dim, 1);
%     for k = 1:I_dim
%         s_k_norm(k) = norm(squeeze(S_full_time(k, k, :)), 2); 
%     end
% 
%     s_k_norm_sq = s_k_norm.^2;
%     VR_k = cumsum(s_k_norm_sq) / sum(s_k_norm_sq);
%     K = find(VR_k >= threshold, 1);
%     if isempty(K), K = I_dim; end
% 
%     % Truncated SVD using calculated K
%     U_fft = zeros(n1, K, n3);   
%     S_fft = zeros(K, K, n3);    
%     V_fft = zeros(n2, K, n3);   
%     for i = 1:half_n3
%         [U_i, S_i, V_i] = svd(X_fft(:,:,i), 'econ');
%         U_fft(:,:,i) = U_i(:, 1:K);
%         S_fft(:,:,i) = S_i(1:K, 1:K);
%         V_fft(:,:,i) = V_i(:, 1:K);
%     end
%     for j = (half_n3 + 1):n3
%         U_fft(:,:,j) = conj(U_fft(:,:, n3 - j + 2));
%         S_fft(:,:,j) = S_fft(:,:, n3 - j + 2);
%         V_fft(:,:,j) = conj(V_fft(:,:, n3 - j + 2));
%     end
% 
%     % Time domain S for clustering
%     S_time = real(ifft(S_fft, [], 3));
% end
% 
% function plot_figure_13_single(s_k_norm, VR_k, patient_id)
%     I_dim = length(s_k_norm);
%     figure('Name', ['Paper Figure 13 - Patient ', num2str(patient_id)], 'Color', 'w', 'Position', [100, 100, 800, 350]);
% 
%     subplot(1, 2, 1);
%     plot(1:I_dim, s_k_norm, 'r--', 'LineWidth', 1.5);
%     xlabel('k', 'FontSize', 12, 'FontAngle', 'italic');
%     ylabel('Magnitude ||s_k||_2', 'FontSize', 12);
%     xlim([0, 1500]);
%     title('Singular Tube Magnitudes');
%     grid on;
% 
%     subplot(1, 2, 2);
%     plot(1:I_dim, VR_k, 'r--', 'LineWidth', 1.5);
%     xlabel('k', 'FontSize', 12, 'FontAngle', 'italic');
%     ylabel('Variance ratio VR(k)', 'FontSize', 12);
%     xlim([0, 1500]); ylim([0, 1.05]);
%     title('Variance Ratio');
%     grid on;
% end
% 
% function idx = cluster_components(S_time, K, num_clusters)
%     [~, ~, n3] = size(S_time);
%     Features = zeros(n3, K); 
%     for k = 1:K
%         Features(:, k) = squeeze(S_time(k,k,:));
%     end
%     % Using K-means per your original script
%     idx = kmeans(Features', num_clusters, 'Distance', 'sqeuclidean', 'Replicates', 10);
% end
% 
% function reconstructed_data = reconstruct_signals(U_fft, S_fft, V_fft, idx, num_clusters, M, N, W, delta)
%     [~, ~, n3] = size(U_fft);
%     reconstructed_data = zeros(num_clusters, M, N);
% 
%     for c = 1:num_clusters
%         component_mask = (idx == c);
%         X_rec_fft = zeros(size(U_fft, 1), size(V_fft, 1), n3);
% 
%         for i = 1:n3
%             S_slice = S_fft(:,:,i);
%             S_slice(~component_mask, ~component_mask) = 0; % Zero out non-cluster values
%             X_rec_fft(:,:,i) = U_fft(:,:,i) * S_slice * V_fft(:,:,i)';
%         end
% 
%         X_rec = real(ifft(X_rec_fft, [], 3));
% 
%         for m = 1:M
%             reconstructed_data(c, m, :) = tblock_diagonal_averaging(X_rec(:,:,m), N, W, delta);
%         end
%     end
% end
% 
% function x_final = block_diagonal_averaging(X_slice, N, W, delta)
%     [~, NumCols] = size(X_slice);
%     x_sum = zeros(1, N);
%     counts = zeros(1, N);
%     for j = 1:NumCols
%         start_idx = (j-1) * delta + 1;
%         end_idx = start_idx + W - 1;
%         if end_idx > N, end_idx = N; end
%         actual_len = end_idx - start_idx + 1;
%         x_sum(start_idx:end_idx) = x_sum(start_idx:end_idx) + X_slice(1:actual_len, j)';
%         counts(start_idx:end_idx) = counts(start_idx:end_idx) + 1;
%     end
%     x_final = x_sum ./ counts;
% end