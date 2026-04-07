function [A, B, Ac, Bc] = initialize_taylor_parameters(para)
    % Initialize Taylor series parameters
    K = para.K;
    K_c = para.K_c;
    
    A = ones(K, K_c);
    B = ones(K, K_c);
    Ac = ones(K, 1);
    Bc = ones(K, 1);
end