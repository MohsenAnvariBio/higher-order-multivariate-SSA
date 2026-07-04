clear; close all; clc;

%% 0. Global Parameters
W = 1500;         
delta = 1;       
num_clusters = 2; 
fs = 200; 
variance_threshold = 0.999; 
channels_to_use = [1, 2, 3, 4, 11, 12]; % M = 6 channels
time_samples = 1:3000;
num_patients = 27;
M = length(channels_to_use);

fprintf('Loading Contaminated_Data.mat into memory...\n');
all_data = load("Contaminated_Data.mat"); 

%% 1. Pre-allocate Storage for Global MAE
% Define the 3 clustering methods we are comparing
cluster_methods = {'M1: K-Means', 'M2: Spectral', 'M3: Hierarchical (Ward''s)'};
num_methods = length(cluster_methods);

% We will store MAE as an array of structs (one struct per method)
empty_struct = struct(...
    'delta', zeros(num_patients, M), ...
    'theta', zeros(num_patients, M), ...
    'alpha', zeros(num_patients, M), ...
    'beta',  zeros(num_patients, M), ...
    'gamma', zeros(num_patients, M));

mae_global = repmat(empty_struct, num_methods, 1);
valid_patients = false(num_patients, 1); % To track which datasets actually exist

%% 2. Batch Processing Loop (All Patients)
fprintf('Starting Batch Evaluation for %d Patients...\n', num_patients);

for p = 1:num_patients
    var_name = sprintf('sim%d_con', p);
    
    if ~isfield(all_data, var_name)
        warning('Variable %s not found. Skipping...', var_name);
        continue;
    end
    
    valid_patients(p) = true;
    
    % Load and normalize data
    raw_data = all_data.(var_name)(channels_to_use, time_samples);
    data = raw_data; % ./ max(abs(raw_data(:))); % Uncomment for global normalization if needed
    [~, N] = size(data);
    
    % Core Pipeline
    X = tbuild_trajectory_tensor(data, W, delta);
    [U_fft, S_fft, V_fft, S_time, K, ~, ~] = tcompute_dynamic_tubal_svd(X, variance_threshold);
    
    % Get all clustering methods
    [idx, idxS, idxC] = tcluster_components(S_time, K, num_clusters);
    
    % Store indices in a cell array for easy looping
    all_idxs = {idx, idxS, idxC};
    
    % --- Loop through each clustering method ---
    for method_idx = 1:num_methods
        current_idx = all_idxs{method_idx};
        
        % Reconstruct signals using the current clustering method
        reconstructed_data_all = treconstruct_signals(U_fft, S_fft, V_fft, current_idx, num_clusters, M, N, W, delta);
        
        % CRITICAL: Identify the Clean Cluster
        clean_cluster_idx = 1; 
        
        % Evaluate MAE for all channels in this patient
        for m = 1:M
            contam_sig = data(m, :);
            recov_sig  = squeeze(reconstructed_data_all(clean_cluster_idx, m, :))';
            
            mae_scores = tevaluate_mae_frequency(contam_sig, recov_sig, fs);
            
            mae_global(method_idx).delta(p, m) = mae_scores.delta;
            mae_global(method_idx).theta(p, m) = mae_scores.theta;
            mae_global(method_idx).alpha(p, m) = mae_scores.alpha;
            mae_global(method_idx).beta(p, m)  = mae_scores.beta;
            mae_global(method_idx).gamma(p, m) = mae_scores.gamma;
        end
    end
    
    fprintf('Processed Patient %d/27\n', p);
end

%% 3. Calculate Global Statistics & Print Table 2 Format
fprintf('\n===============================================================================================================\n');
fprintf('Table 2 Recreation: Global MAE Results (Mean +/- Std) over %d patients and %d channels\n', sum(valid_patients), M);
fprintf('===============================================================================================================\n');
fprintf('%-18s | %-14s | %-14s | %-14s | %-14s | %-14s\n', 'Clustering Method', 'delta', 'theta', 'alpha', 'beta', 'low-gamma');
fprintf('---------------------------------------------------------------------------------------------------------------\n');

bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};

% Loop through methods and print a row for each
for method_idx = 1:num_methods
    results_str = cell(1, 5);
    
    for b = 1:length(bands)
        band_name = bands{b};
        
        % Get valid data for the current method and flatten to 1D
        valid_data_matrix = mae_global(method_idx).(band_name)(valid_patients, :);
        flat_data = valid_data_matrix(:); 
        
        avg_val = mean(flat_data, 'omitnan');
        std_val = std(flat_data, 'omitnan');
        
        results_str{b} = sprintf('%6.3f ± %5.3f', avg_val, std_val);
    end
    
    % Print the formatted row for the specific method
    fprintf('%-18s | %-14s | %-14s | %-14s | %-14s | %-14s\n', ...
        cluster_methods{method_idx}, results_str{1}, results_str{2}, results_str{3}, results_str{4}, results_str{5});
end
fprintf('===============================================================================================================\n\n');

%% 4. Dynamic Visualization
[my_means, my_stds, comparison_fig] = tcomparison_mae(mae_global, valid_patients, cluster_methods);