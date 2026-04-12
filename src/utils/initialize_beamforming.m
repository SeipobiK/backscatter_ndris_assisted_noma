function W_init = initialize_beamforming(para, H_eff)
    % Initialize beamforming vectors using maximum ratio transmission
    K = para.K;
    M = para.M;
    
    W_init = zeros(M, K);
    
    for k = 1:K
        h = H_eff{k, 2};  % Use primary subcarrier
        if norm(h) > 0
            W_init(:, k) = h / norm(h);
        else
            W_init(:, k) = randn(M, 1) + 1i*randn(M, 1);
            W_init(:, k) = W_init(:, k) / norm(W_init(:, k));
        end
    end
end