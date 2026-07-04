% =========================================================================
% 1. User Input & Parameter Initialization
% =========================================================================
clc, clear, close all

patient_id = input('Enter the dataset number to process (1 to 27): ');

% Input validation
if isempty(patient_id) || patient_id < 1 || patient_id > 27
    fprintf('Invalid input. Defaulting to dataset 3.\n');
    patient_id = 3; 
end

fs = 200; 
channels_to_use = [1, 2, 3, 4, 11, 12]; % Fp1, Fp2, F3, F4, F7, F8
time_samples = 1:5800;

all_data = load("Contaminated_Data.mat"); 
var_name = sprintf('sim%d_con', patient_id);

if ~isfield(all_data, var_name)
    error('Variable %s not found in Contaminated_Data.mat.', var_name);
end

raw_data = all_data.(var_name)(channels_to_use, time_samples);
[M, N] = size(raw_data); 

% CRITICAL FIX 1: Normalize globally, NOT channel-by-channel.
data = raw_data; % ./ max(abs(raw_data(:)));

% Calculate the time vector in seconds based on your time_samples and fs
time = (0:length(time_samples)-1) / fs; 

channel_names = {'Fp1', 'Fp2', 'F3', 'F4', 'F7', 'F8'};
num_channels = length(channel_names);

% =========================================================================
% 3. Plotting (Figure 12 Style)
% =========================================================================
figure('Name','EG Signals Mixed with EOG Artifacts' , 'Color', 'w', 'Position', [100, 100, 900, 400]);

% Dynamically update the title with the requested patient ID
sgtitle(['EEG Signals Mixed with EOG Artifacts: Patient ', num2str(patient_id)], ...
    'FontSize', 16, 'FontWeight', 'bold');

for i = 1:num_channels
    subplot(num_channels, 1, i);
    plot(time, raw_data(i, :), 'k', 'LineWidth', 0.8); 
    
    xlim([time(1) time(end)]);
    
    % Formatting: clean "floating line" aesthetic
    box off;
    set(gca, 'YTick', []); 
    set(gca, 'YColor', 'none'); % This hides the spine (and the ylabel!)
    
    % FIX: Use text() to manually place the channel name on the left
    y_limits = ylim; % Get current y-axis boundaries to find the vertical center
    x_offset = time(1) - 0.02 * (time(end) - time(1)); % Push text slightly left of x=0
    
    text(x_offset, mean(y_limits), channel_names{i}, ...
        'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    
    if i < num_channels
        % Hide x-axis for all upper subplots
        set(gca, 'XTick', []);
        set(gca, 'XColor', 'none'); 
    else
        % Show x-axis only on the very bottom plot
        set(gca, 'XColor', 'k'); 
        xlabel('Time (seconds)', 'FontSize', 12, 'FontWeight', 'bold');
    end
end