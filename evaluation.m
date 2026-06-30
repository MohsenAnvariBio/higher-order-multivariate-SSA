%% 6. Evaluation: MAE and Correlation in the Frequency Domain
% Note: Since HO-MSSA separates the signal into clusters, you must specify 
% which cluster represents the "recovered/clean" EEG. Update this variable
% based on your visual inspection of the plots.
target_cluster = 1; 

% Frequency bands defined in the paper (excluding Delta 0.5-4 Hz)
bands = struct();
bands.Delta = [0.5, 4];
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

% Correlation window parameters: 2*delta_f = 3 Hz -> delta_f = 1.5 Hz
delta_f = 1.5; 
delta_bins = round(delta_f / df);   % Convert Hz to number of FFT bins

disp('Calculating Frequency Domain Metrics...');

for m = 1:M
    % 1. Extract Contaminated (x) and Recovered (y) signals
    x = data(m, :); 
    y = squeeze(reconstructed_data_all(target_cluster, m, :))';
    
    % 2. Calculate Fourier Coefficients
    X_k = fft(x);
    Y_k = fft(y);
    
    % Extract only the positive frequencies to match the paper's spectrum
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
    
    % % --- Correlation rho(f) Calculation (Equation 22) ---
    % rho_f = zeros(1, half_idx);
    % for i = 1:half_idx
    %     % Define the moving window boundaries [f - delta_f, f + delta_f]
    %     start_idx = max(1, i - delta_bins);
    %     end_idx = min(half_idx, i + delta_bins);
    % 
    %     X_win = X_k_half(start_idx:end_idx);
    %     Y_win = Y_k_half(start_idx:end_idx);
    % 
    %     % Numerator terms
    %     term1 = conj(X_win) .* Y_win;
    %     term2 = X_win .* conj(Y_win);
    %     numerator = 0.5 * abs(sum(term1 + term2));
    % 
    %     % Denominator terms
    %     sum_xx = sum(X_win .* conj(X_win));
    %     sum_yy = sum(Y_win .* conj(Y_win));
    %     denominator = sqrt(abs(sum_xx * sum_yy));
    % 
    %     % Avoid division by zero
    %     if denominator == 0
    %         rho_f(i) = 0;
    %     else
    %         rho_f(i) = numerator / denominator;
    %     end
    % end
    % 
    % % --- Plotting Correlation rho(f) ---
    % % Plotting just for the first channel as an example, matching Figure 14(e/f)
    % if m == 1 
    %     figure('Name', 'Frequency Domain Correlation \rho(f)', 'Color', 'w');
    %     plot(f_half, rho_f, 'b', 'LineWidth', 1.5);
    %     xlim([0 40]); % The paper explicitly evaluates the spectrum up to 40 Hz
    %     ylim([0 1.1]);
    %     xlabel('Frequency (Hz)');
    %     ylabel('\rho(f)');
    %     title(['Correlation \rho(f) for Channel ', num2str(m)]);
    %     grid on;
    % end
end

% --- Display the MAE Results Table in the Command Window ---
disp('--- Mean Absolute Error (MAE) Results ---');
RowNames = arrayfun(@(x) sprintf('Channel_%d', x), 1:M, 'UniformOutput', false);
MAE_Table = array2table(MAE_results, 'VariableNames', band_names, 'RowNames', RowNames);
disp(MAE_Table);