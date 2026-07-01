function [beta, zeta, y, z, sum_rate_const_alpha, sum_rate_const_eta, duals_vars, duals_vars_c] = ...
    compute_fp_gual_vars(para, H, H_c, decoding_order, alpha, w_k, eta, ...
    duals_vars_init, learning_rate, duals_vars_init_c)

    %% ---- Parameters ----
    K     = para.K;
    K_c   = para.K_c;
    noise = para.noise;

    %% ---- Initialise ----
    A     = zeros(K, K_c);
    B     = zeros(K, K_c);
  duals_vars   = duals_vars_init;    % ← start from input, not zero
duals_vars_c = duals_vars_init_c;  % ← same
    A_c   = zeros(K, 1);
    B_c   = zeros(K, 1);

    beta  = zeros(K, K_c);
    I_tot = zeros(K, K_c);
    y     = zeros(K, K_c);
    zeta  = zeros(K, 1);
    z     = zeros(K, 1);
    
    gamma_noma = 2^para.R_min_n - 1;
    gamma_back = 2^para.R_min_c - 1;

    sum_rate_const_alpha = 0;
    sum_rate_const_eta = 0;
    
    %% ===== MAIN LOOP =====
    for k = 1:K
        order_k = decoding_order(k, :);
        strong_user = order_k(end);

        for i = 1:K_c
            %% Interference calculations (your existing code)
            inter = 0;
            inter_b = 0;
            
            for j = 1:K
                if j ~= k
                    inter = inter + abs(H{k,i} * w_k(:,j)).^2;
                end
                
                inter_bb = inter_b + abs(H_c{k,i} * w_k(:,j)).^2 * eta(j);
            end

            intra = 0;
            pos = find(order_k == i);
            H_ki = abs(H{k,i} * w_k(:,k)).^2;

            for idx = pos+1:K_c
                j_user = order_k(idx);
                intra = intra + alpha(k, j_user) * H_ki;
            end

            A(k,i) = H_ki * alpha(k,i);
            B(k,i) = intra + inter + inter_b + noise;
            
            % Dual variable update (fixed typo)
            duals_vars(k,i) = max((duals_vars_init(k,i) + learning_rate * (A(k,i) - gamma_noma * B(k,i))), 0);

            I_tot(k,i) = inter + noise;
            beta(k,i) = A(k,i) / B(k,i);
            sum_rate_const_alpha = sum_rate_const_alpha + log2(1+beta(k,i)) - beta(k,i);
            y(k,i) = sqrt((1 + beta(k,i)) * A(k,i)) / (A(k,i) + B(k,i));

            %% Backscatter for strong user
            if i == strong_user
                H_c_kk = abs(H_c{k,i} * w_k(:,k)).^2;
                A_c(k) = H_c_kk * eta(k);
                inter_b_others = inter_b - H_c_kk * eta(k);
                B_c(k) = inter_b_others + inter + noise;
                
                duals_vars_c(k) = max((duals_vars_init_c(k) + learning_rate * (A_c(k) - gamma_back * B_c(k))), 0);

                zeta(k) = A_c(k) / B_c(k);
                sum_rate_const_eta = sum_rate_const_eta + log2(1+zeta(k)) - zeta(k);
                z(k) = sqrt((1 + zeta(k)) * A_c(k)) / (A_c(k) + B_c(k));
            end
        end
    end
end