clear; close all; clc;

%% 1. Define the Data (Correctly Aligned: Delta, Theta, Alpha, Beta, Gamma)
bands = {'Delta (0.5-4 Hz)', 'Theta (4-8 Hz)', 'Alpha (8-12 Hz)', 'Beta (12-30 Hz)', 'Low-Gamma (30-40 Hz)'};
methods = {'Paper (HO-MSSA)', 'Standard (idx)', 'Symmetric (idxS)', 'Custom (idxC)'};

% Means Matrix: [Methods x Bands]
means = [
    102.1,  1.301,  0.705,  0.133,  0.051;  % Paper (HO-MSSA)
    37.267, 0.261,  3.499,  0.044,  0.034;  % Standard (idx)
    37.302, 2.826,  9.174,  0.952,  0.399;  % Symmetric (idxS)
    38.152, 2.744,  10.598, 0.970,  0.419   % Custom (idxC)
    ];

% Standard Deviations Matrix: [Methods x Bands]
stds = [
    37.78,  0.374,  0.462,  0.085,  0.029;  % Paper (HO-MSSA)
    36.020, 1.044,  5.470,  0.066,  0.015;  % Standard (idx)
    34.948, 1.716,  7.644,  0.773,  0.421;  % Symmetric (idxS)
    36.180, 2.779,  10.508, 1.300,  0.581   % Custom (idxC)
    ];

%% 2. Plotting: Subplots per Frequency Band
figure('Name', 'MAE Comparison by Frequency Band', 'Position', [100, 100, 1200, 400]);
colors = [0.2 0.5 0.7;  % Blue (Paper)
    0.3 0.7 0.4;  % Green (Standard)
    0.9 0.6 0.2;  % Orange (Symmetric)
    0.8 0.3 0.3]; % Red (Custom)

for b = 1:5
    subplot(1, 5, b);

    % Get data for this specific band
    band_means = means(:, b);
    band_stds = stds(:, b);

    % Create grouped bars for the 4 methods
    for m = 1:4
        bar(m, band_means(m), 'FaceColor', colors(m,:));
        hold on;
    end

    % Add error bars
    errorbar(1:4, band_means, band_stds, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

    % Formatting
    title(bands{b}, 'FontSize', 12, 'FontWeight', 'bold');
    xticks(1:4);
    xticklabels({'Paper', 'Standard', 'Sym', 'Custom'});
    xtickangle(45);
    grid on;

    % Only label the Y-axis on the far left plot to keep it clean
    if b == 1
        ylabel('Mean Absolute Error (MAE)', 'FontSize', 12, 'FontWeight', 'bold');
    end
end

% Add a global legend to the first subplot
subplot(1, 5, 1);
legend(methods, 'Location', 'northeast', 'FontSize', 10);