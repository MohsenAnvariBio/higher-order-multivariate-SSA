function [my_means, my_stds, hFig] = tcomparison_mae(mae_global, valid_patients, cluster_methods)
% TPLOT_DYNAMIC_MAE_COMPARISON Visualizes your MAE results vs the paper benchmark.
%
% Inputs:
%   mae_global      - Struct array containing the MAE results for each method.
%   valid_patients  - Logical array indicating which patients were processed.
%   cluster_methods - Cell array of strings with the names of your methods.
%
% Outputs:
%   my_means        - Matrix of calculated means [Methods x Bands]
%   my_stds         - Matrix of calculated standard deviations [Methods x Bands]
%   hFig            - Handle to the generated figure

% 1. Define the fixed benchmark data from the paper (HO-MSSA)
% Order: [Delta, Theta, Alpha, Beta, Low-Gamma]
paper_means = [102.1,  1.301,  0.705,  0.133,  0.051];
paper_stds  = [37.78,  0.374,  0.462,  0.085,  0.029];

bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
num_methods = length(cluster_methods);

% 2. Dynamically extract your calculated means and stds into matrices
my_means = zeros(num_methods, 5);
my_stds  = zeros(num_methods, 5);

for m_idx = 1:num_methods
    for b_idx = 1:length(bands)
        band_name = bands{b_idx};

        % Get valid data and flatten to 1D
        valid_data_matrix = mae_global(m_idx).(band_name)(valid_patients, :);
        flat_data = valid_data_matrix(:); 

        % Compute dynamic statistics
        my_means(m_idx, b_idx) = mean(flat_data, 'omitnan');
        my_stds(m_idx, b_idx)  = std(flat_data, 'omitnan');
    end
end

% 3. Combine Paper data (Row 1) with Your data (Rows 2+)
all_means = [paper_means; my_means];
all_stds  = [paper_stds; my_stds];

% 4. Plotting Setup
plot_bands = {'Delta (0.5-4 Hz)', 'Theta (4-8 Hz)', 'Alpha (8-12 Hz)', 'Beta (12-30 Hz)', 'Low-Gamma (30-40 Hz)'};
plot_methods = [{'Paper (HO-MSSA)'}, cluster_methods]; % Prepend paper to legend

hFig = figure('Name', 'Dynamic MAE Comparison', 'Position', [100, 100, 1200, 450]);
sgtitle('Performance Comparison: EOG Artifact Removal across Frequency Bands', 'FontSize', 16, 'FontWeight', 'bold');

% Define colors: Paper is Blue, followed by Green, Orange, Red, etc.
colors = [0.2 0.5 0.7;  
    0.3 0.7 0.4;  
    0.9 0.6 0.2;  
    0.8 0.3 0.3;
    0.5 0.3 0.7]; 

% Ensure we don't run out of colors if you add more methods
colors = repmat(colors, ceil((num_methods + 1) / size(colors, 1)), 1);

for b = 1:5
    subplot(1, 5, b);

    band_means = all_means(:, b);
    band_stds = all_stds(:, b);

    % Create grouped bars
    for m = 1:(num_methods + 1)
        bar(m, band_means(m), 'FaceColor', colors(m,:));
        hold on;
    end

    % Add error bars
    errorbar(1:(num_methods + 1), band_means, band_stds, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

    % Formatting
    title(plot_bands{b}, 'FontSize', 12, 'FontWeight', 'bold');
    xticks(1:(num_methods + 1));
    xticklabels({'Paper', 'M1', 'M2', 'M3', 'M4', 'M5'}); % Shortened for spacing
    xtickangle(45);
    grid on;

    if b == 1
        ylabel('Mean Absolute Error (MAE)', 'FontSize', 12, 'FontWeight', 'bold');
    end
end

% Add a global legend to the first subplot
subplot(1, 5, 1);
legend(plot_methods, 'Location', 'best', 'FontSize', 10);

end