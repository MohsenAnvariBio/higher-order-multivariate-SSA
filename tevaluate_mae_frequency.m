function mae_results = tevaluate_mae_frequency(x_con, x_rec, fs)
bands = struct('delta', [0.5, 4], 'theta', [4, 8], 'alpha', [8, 12], 'beta', [12, 30], 'gamma', [30, 40]);
band_names = fieldnames(bands);

window_length = min(length(x_con), 2 * fs); 
[P_con, f] = pwelch(x_con, window_length, [], [], fs);
[P_rec, ~] = pwelch(x_rec, window_length, [], [], fs);

mae_results = struct();
for i = 1:length(band_names)
    b_name = band_names{i};
    b_range = bands.(b_name);
    idx = find(f >= b_range(1) & f <= b_range(2));

    if isempty(idx)
        mae_results.(b_name) = NaN; 
        continue;
    end

    abs_diff = abs(P_con(idx) - P_rec(idx));
    mae_results.(b_name) = sum(abs_diff) / length(idx);
end
end