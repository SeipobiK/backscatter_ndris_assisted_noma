function W = extract_beamforming_vectors(W_opt)
    % Extract beamforming vectors from optimization results
    K = size(W_opt, 3);
    M = size(W_opt, 1);
    W = zeros(M, K);
    
    for k = 1:K
        % Assuming max_eigVect is a custom function that returns [eigenvector, eigenvalue]
        [W_max, max_eigenvalue] = max_eigVect(W_opt(:, :, k));
        W(:, k) = sqrt(max_eigenvalue) * W_max;
        % disp(['Max eigenvalue for user ', num2str(k), ': ', num2str(max_eigenvalue)]);
    end
end