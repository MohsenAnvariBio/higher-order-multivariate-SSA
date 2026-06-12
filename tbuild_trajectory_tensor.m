function X = tbuild_trajectory_tensor(data, W, delta)
[M, N] = size(data);
num_cols = floor((N - W) / delta) + 1;
X = zeros(W, num_cols, M);
for m = 1:M
    for j = 1:num_cols
        start_idx = (j-1)*delta + 1;
        X(:, j, m) = data(m, start_idx : start_idx + W - 1);
    end
end
end