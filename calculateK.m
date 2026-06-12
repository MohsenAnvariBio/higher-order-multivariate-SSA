clear; close all; clc;

%% 0. Global Parameters
W = 1500;         
delta = 1;       
fs = 200; 
num_patients = 27; 
variance_threshold = 0.999; % 99.9% for calculating K
channels_to_use = [1, 2, 3, 4, 11, 12]; % Fp1, Fp2, F3, F4, F7, F8
time_samples = 1:3000;

% --- DATA LOADING ---
% Load the entire file into a structure variable
fprintf('Loading Contaminated_Data.mat into memory...\n');
all_data = load("Contaminated_Data.mat"); 

%% 1. Process 27 Patients
% We initialize the storage arrays based on the dimensions we expect.
% Number of columns: floor((3000 - 1500) / 1) + 1 = 1501
I_dim = min(W, floor((length(time_samples) - W) / delta) + 1);
all_s_k_norms = zeros(num_patients, I_dim);
all_VR_k = zeros(num_patients, I_dim);
all_K = zeros(num_patients, 1);

fprintf('Processing %d patients...\n', num_patients);

for p = 1:num_patients
    % 1a. Dynamically construct the variable name (e.g., 'sim1_con', 'sim2_con')
    var_name = sprintf('sim%d_con', p);

    % 1b. Check if the variable actually exists in the loaded file
    if ~isfield(all_data, var_name)
        warning('Variable %s not found in Contaminated_Data.mat. Skipping...', var_name);
        continue;
    end

    % 1c. Extract the specific channels and time samples for this patient
    patient_data = all_data.(var_name)(channels_to_use, time_samples);

    % CRITICAL FIX: Normalize globally for each patient
    patient_data = patient_data ./ max(abs(patient_data(:)));

    % --- STEP 1: Construct Trajectory Tensor ---
    X = tbuild_trajectory_tensor(patient_data, W, delta);

    % --- STEP 2: Compute Tubal Norms, VR(k), and K ---
    [s_k_norm, VR_k, K] = tcompute_tubal_metrics(X, variance_threshold);

    % Store the metrics for plotting later
    all_s_k_norms(p, :) = s_k_norm;
    all_VR_k(p, :) = VR_k;
    all_K(p) = K;

    fprintf('Patient %2d (%s): Calculated K = %d\n', p, var_name, K);
end

%% 2. Display Average K Statistics (Matches Table 1 logic in paper)
% Filter out any empty rows in case some variables were missing and skipped
valid_idx = all_K > 0; 
valid_K = all_K(valid_idx);

mean_K = mean(valid_K);
std_K = std(valid_K);
fprintf('\n----------------------------------------\n');
fprintf('Tubal Rank K over %d valid patients: %.2f ± %.2f\n', sum(valid_idx), mean_K, std_K);
fprintf('----------------------------------------\n');

%% 3. Plot Figure 13 (Averaged over 27 patients)
% Only plot using the rows where data was successfully processed
tplot_figure_13(all_s_k_norms(valid_idx, :), all_VR_k(valid_idx, :), I_dim);
