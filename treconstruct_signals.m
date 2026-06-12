function reconstructed_data = treconstruct_signals(U_fft, S_fft, V_fft, idx, num_clusters, M, N, W, delta)
[~, ~, n3] = size(U_fft);
reconstructed_data = zeros(num_clusters, M, N);

for c = 1:num_clusters
    component_mask = (idx == c);
    X_rec_fft = zeros(size(U_fft, 1), size(V_fft, 1), n3);

    for i = 1:n3
        S_slice = S_fft(:,:,i);
        S_slice(~component_mask, ~component_mask) = 0; % Zero out non-cluster values
        X_rec_fft(:,:,i) = U_fft(:,:,i) * S_slice * V_fft(:,:,i)';
    end

    X_rec = real(ifft(X_rec_fft, [], 3));

    for m = 1:M
        reconstructed_data(c, m, :) = tblock_diagonal_averaging(X_rec(:,:,m), N, W, delta);
    end
end
end