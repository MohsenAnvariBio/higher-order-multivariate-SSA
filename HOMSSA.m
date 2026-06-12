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
