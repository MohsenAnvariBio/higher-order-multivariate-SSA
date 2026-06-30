%% 6. Evaluation: MAE and Correlation in the Frequency Domain
% Note: Since HO-MSSA separates the signal into clusters, you must specify 
% which cluster represents the "recovered/clean" EEG. Update this variable
% based on your visual inspection of the plots.
target_cluster = 1; 

% Frequency bands defined in the paper 
bands = struct();
bands.Delta    = [0.5, 4];
bands.Theta    = [4, 8];
bands.Alpha    = [8, 12];
bands.Beta     = [12, 30];
bands.LowGamma = [30, 40];
band_names = fieldnames(bands);

% Frequency axis setup
df = fs / N;                        % Frequency resolution (Hz per bin)
f_axis = (0:N-1) * df;              % Full frequency axis
half_idx = floor(N/2);              % Nyquist limit index
f_half = f_axis(1:half_idx);        % Positive frequencies

% Pre-allocate MAE results matrix
MAE_results = zeros(M, length(band_names));

disp('Calculating Frequency Domain Metrics...');

for m = 1:M
    % 1. Extract Contaminated (x) and Recovered (y) signals
    x = data(m, :); 
    y = squeeze(reconstructed_data_all(target_cluster, m, :))';
    
    % 2. Calculate Fourier Coefficients
    X_k = fft(x);
    Y_k = fft(y);
    
    % Extract only the positive frequencies
    X_k_half = X_k(1:half_idx);
    Y_k_half = Y_k(1:half_idx);
    
    % 3. Calculate Power Spectral Density (PSD)
    P_x = (1/(fs*N)) * abs(X_k_half).^2;
    P_y = (1/(fs*N)) * abs(Y_k_half).^2;
    
    % --- MAE Calculation (Equation 21) ---
    for b = 1:length(band_names)
        band_range = bands.(band_names{b});
        % Find frequency bins that fall inside the current band
        band_idx = find(f_half >= band_range(1) & f_half <= band_range(2));
        
        if ~isempty(band_idx)
            % MAE = Sum(|P_con - P_rec|) / |B_f|
            mae_val = sum(abs(P_x(band_idx) - P_y(band_idx))) / length(band_idx);
            MAE_results(m, b) = mae_val;
        end
    end
    
    % --- Plotting PSD and Frequency Bands (For Channel 1) ---
    % Change "m == 1" to "true" if you want to generate a plot for every channel
    if m == 6 
        figure('Name', ['PSD and Frequency Bands: Channel ', num2str(m)], 'Color', 'w');
        hold on;
        
        % Define colors for the shaded frequency bands
        band_colors = lines(length(band_names));
        
        % Calculate a dynamic Y-limit based on the maximum power up to 45 Hz
        plot_idx = f_half <= 45; 
        y_max = max(max(P_x(plot_idx)), max(P_y(plot_idx))) * 1.1; 
        if y_max == 0, y_max = 1; end % Prevent axis error if signals are flat
        
        % Draw shaded rectangles for each frequency band
        for b = 1:length(band_names)
            band_range = bands.(band_names{b});
            patch([band_range(1) band_range(2) band_range(2) band_range(1)], ...
                  [0 0 y_max y_max], band_colors(b,:), ...
                  'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                  
            % Add band labels at the top of the shaded regions
            mid_f = mean(band_range);
            text(mid_f, y_max * 0.95, band_names{b}, 'HorizontalAlignment', 'center', ...
                 'FontSize', 9, 'FontWeight', 'bold', 'Color', band_colors(b,:)*0.7);
        end
        
        % Plot the actual PSD lines
        h1 = plot(f_half, P_x, 'r', 'LineWidth', 1.2);
        h2 = plot(f_half, P_y, 'b', 'LineWidth', 1.2);
        
        % Formatting
        xlim([0 45]); % View slightly past the 40 Hz Low-Gamma cutoff
        ylim([0 y_max]);
        xlabel('Frequency (Hz)');
        ylabel('Power Spectral Density (PSD)');
        title(['Frequency Domain Separation - Channel ', num2str(m)]);
        legend([h1, h2], {'Contaminated (P_{con})', 'Recovered (P_{rec})'}, 'Location', 'northeast');
        grid on;
        hold off;
    end
end

% --- Display the MAE Results Table in the Command Window ---
disp('--- Mean Absolute Error (MAE) Results ---');
RowNames = arrayfun(@(x) sprintf('Channel_%d', x), 1:M, 'UniformOutput', false);
MAE_Table = array2table(MAE_results, 'VariableNames', band_names, 'RowNames', RowNames);
disp(MAE_Table);