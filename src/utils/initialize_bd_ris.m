function Theta = initialize_bd_ris(para, type)
    % INITIALIZE_BD_RIS - Initialize fully-connected BD-RIS scattering matrix
    %
    % Inputs:
    %   para  - structure with .N (number of RIS elements)
    %   type  - initialization type: 'eye', 'random', 'hadamard', 'fourier'
    %
    % Output:
    %   Theta - N x N unitary matrix (fully-connected BD-RIS)
    
    N = para.N;
    
    switch type
        case 'eye'
            % Identity matrix (no reflection/transmission)
            Theta = eye(N);
            
        case 'random'
            % Random unitary matrix via QR decomposition
            Z = (randn(N, N) + 1i*randn(N, N)) / sqrt(2);
            [Q, ~] = qr(Z);
            Theta = Q;
            
        case 'hadamard'
            % Hadamard-based unitary (real, good for large N)
            H = hadamard(N);
            Theta = H / sqrt(N);
            
        case 'fourier'
            % Discrete Fourier Transform matrix
            F = dftmtx(N);
            Theta = F / sqrt(N);
            
        otherwise
            % Default: identity
            Theta = eye(N);
    end

    % Verify unitary property
    if norm(Theta * Theta' - eye(N), 'fro') > 1e-10
        warning('Initial Theta is not unitary. Normalizing...');
        [U, ~, V] = svd(Theta);
        Theta = U * V';
    end
end