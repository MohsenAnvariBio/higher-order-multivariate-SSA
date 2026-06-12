
function x_final = tblock_diagonal_averaging(X_slice, N, W, delta)
[~, NumCols] = size(X_slice);
x_sum = zeros(1, N);
counts = zeros(1, N);
for j = 1:NumCols
    start_idx = (j-1) * delta + 1;
    end_idx = start_idx + W - 1;
    if end_idx > N, end_idx = N; end
    actual_len = end_idx - start_idx + 1;
    x_sum(start_idx:end_idx) = x_sum(start_idx:end_idx) + X_slice(1:actual_len, j)';
    counts(start_idx:end_idx) = counts(start_idx:end_idx) + 1;
end
x_final = x_sum ./ counts;
end