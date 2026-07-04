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
% We will store MAE as [Patients x Channels] for each band
mae_global = struct(...
    'delta', zeros(num_patients, M), ...
    'theta', zeros(num_patients, M), ...
    'alpha', zeros(num_patients, M), ...
    'beta',  zeros(num_patients, M), ...
    'gamma', zeros(num_patients, M));
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
    data = raw_data ;%./ max(abs(raw_data(:))); % Global normalization
    [~, N] = size(data);
    
    % Pipeline
    X = tbuild_trajectory_tensor(data, W, delta);
    [U_fft, S_fft, V_fft, S_time, K, ~, ~] = tcompute_dynamic_tubal_svd(X, variance_threshold);
    idx = tcluster_components(S_time, K, num_clusters);
    reconstructed_data_all = treconstruct_signals(U_fft, S_fft, V_fft, idx, num_clusters, M, N, W, delta);
    
    % --- CRITICAL: Identify the Clean Cluster ---
    % Defaulting to 1 for this example. If you have an automated way to 
    % detect the artifact cluster (e.g., max variance in 0-4Hz), apply it here.
    clean_cluster_idx = 1; 
    
    % Evaluate MAE for all channels in this patient
    for m = 1:M
        contam_sig = data(m, :);
        recov_sig  = squeeze(reconstructed_data_all(clean_cluster_idx, m, :))';
        
        mae_scores = tevaluate_mae_frequency(contam_sig, recov_sig, fs);
        
        mae_global.delta(p, m) = mae_scores.delta;
        mae_global.theta(p, m) = mae_scores.theta;
        mae_global.alpha(p, m) = mae_scores.alpha;
        mae_global.beta(p, m)  = mae_scores.beta;
        mae_global.gamma(p, m) = mae_scores.gamma;
    end
    
    fprintf('Processed Patient %d/27\n', p);
end

%% 3. Calculate Global Statistics & Print Table 2 Format
fprintf('\n========================================================================================\n');
fprintf('Table 2 Recreation: Global MAE Results (Mean +/- Std) over %d patients and %d channels\n', sum(valid_patients), M);
fprintf('========================================================================================\n');
fprintf('%-12s | %-12s | %-12s | %-12s | %-12s\n', 'delta', 'theta', 'alpha', 'beta', 'low-gamma');
fprintf('----------------------------------------------------------------------------------------\n');

% Extract the valid data, flatten the [Patients x Channels] matrix into a 1D array, 
% and calculate Mean +/- Std.
bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
results_str = cell(1, 5);
for b = 1:length(bands)
    band_name = bands{b};
    
    % Get data only for patients that were successfully loaded
    valid_data_matrix = mae_global.(band_name)(valid_patients, :);
    
    % Flatten matrix to a single column vector representing ALL channels globally
    flat_data = valid_data_matrix(:); 
    
    avg_val = mean(flat_data, 'omitnan');
    std_val = std(flat_data, 'omitnan');
    
    results_str{b} = sprintf('%6.3f ± %5.3f', avg_val, std_val);
end

% Print the formatted row
fprintf('%-12s | %-12s | %-12s | %-12s | %-12s\n', ...
    results_str{1}, results_str{2}, results_str{3}, results_str{4}, results_str{5});
fprintf('========================================================================================\n\n');

%% 4. Plot Power Spectral Density (PSD) matching Reference Image
fprintf('Generating PSD plot for the last processed signal...\n');

% Compute PSD using Welch's method for both signals
[pxx_c, f_psd] = pwelch(contam_sig, [], [], [], fs);
[pxx_r, ~]     = pwelch(recov_sig, [], [], [], fs);

% Convert to Decibels (dB/Hz)
psd_c_dB = 10*log10(pxx_c);
psd_r_dB = 10*log10(pxx_r);

% Define Frequency Bands: [Start Hz, End Hz]
bands_limits = [
    0, 4;   % Delta
    4, 8;   % Theta
    8, 13;  % Alpha
    13, 30; % Beta
    30, 100 % Gamma
];

% Colors for frequency bands
band_colors = [
    0, 0.4470, 0.7410;       % Blue
    0.8500, 0.3250, 0.0980;  % Orange
    0.9290, 0.6940, 0.1250;  % Yellow
    0.4940, 0.1840, 0.5560;  % Purple
    0.4660, 0.6740, 0.1880   % Green
];

figure('Name', 'Power Spectral Density', 'Position', [100 100 900 450]);
hold on;

% Determine Y-axis limits dynamically based on data
ymin = min([psd_c_dB; psd_r_dB]) - 5;
ymax = max([psd_c_dB; psd_r_dB]) + 5;

% 1. Draw the colored bands in the background first
for b = 1:size(bands_limits, 1)
    x_patch = [bands_limits(b,1) bands_limits(b,2) bands_limits(b,2) bands_limits(b,1)];
    y_patch = [ymin ymin ymax ymax];
    
    % Draw transparent patch (HandleVisibility off so it doesn't crowd the legend)
    patch(x_patch, y_patch, band_colors(b,:), 'FaceAlpha', 0.15, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% 2. Plot the PSD Lines on top (Black for Contaminated, Blue for Denoised)
plot(f_psd, psd_c_dB, 'k', 'LineWidth', 1.2, 'DisplayName', 'Contaminated EEG');
plot(f_psd, psd_r_dB, 'b', 'LineWidth', 1.2, 'DisplayName', 'Denoised EEG');
title('Power spectral density', 'FontWeight', 'normal', 'FontSize', 16);
% 3. Formatting to match the image
xlabel('Frequency (Hz)', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Power/Frequency (dB/Hz)', 'FontWeight', 'bold', 'FontSize', 16);
grid on;
xlim([0 100]); % Cap X-axis at 100Hz as shown in your image
ylim([ymin ymax]); 
legend('Location', 'northeast', 'FontSize', 16);

hold off;