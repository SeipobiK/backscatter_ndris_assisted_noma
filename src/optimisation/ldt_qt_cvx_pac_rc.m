function [p_alloc, eta, obj] = ldt_qt_cvx_pac_rc(...
    para, H, H_c,h_c_eh, I_tot, beta, t, zeta, s, decoding_order, p_init, eta_init)

    % Parameters
    K = para.K;
    K_c = para.K_c;
    R_min_noma = 2^para.R_min_n - 1;      % SINR target for NOMA users
    R_min_back = 2^para.R_min_c - 1;      % SINR target for backscatter
    eh=para.bst_threshold;
    rho=para.rho;
    
    % Initialize with previous values if provided
    if nargin < 10
        p_init = ones(K, K_c) / K_c;
    end
    if nargin < 11
        eta_init = 0.5 * ones(K, 1);
    end

    cvx_begin quiet
        % cvx_solver mosek
        
        % Variables (using p for power allocation to avoid conflicts)
        variable p(K, K_c) nonnegative
        variable eta(K, 1) nonnegative
        variable delta_G nonnegative
        
        % Objective
        sum_rate = 0;
        
        for k = 1:K
            order_k = decoding_order(k, :);  % weak → strong
            strong_user = order_k(end);       % strongest user (decodes last)
            
            for i = 1:K_c
                % Inter-backscatter interference from other clusters
                inter_b = sum(eta) - eta(k);  % Sum of all eta except own cluster
                
                % Intra-cluster interference (from stronger users)
                pos = find(order_k == i);
                intra = 0;
                for idx = pos+1:K_c
                    j_user = order_k(idx);
                    intra = intra + p(k, j_user);
                end
                
                % Objective: NOMA user rate (using LDT + QT)
                % sum_rate = sum_rate + log2(1 + beta(k,i)) - beta(k,i) ...
                %     + 2 * t(k,i) * sqrt((1 + beta(k,i)) * H(k,i) * p(k,i)) ...
                %     - t(k,i)^2 * (I_tot(k,i) + H_c(k,i) * (inter_b + eta(k)) + H(k,i) * (intra + p(k,i)));
                
                % QoS constraint for NOMA user
                H(k,i) * p(k,i)+delta_G >= R_min_noma * (I_tot(k,i) + H_c(k,i) * (inter_b + eta(k)) + H(k,i) * intra);
                
                % Backscatter link (only for strongest user)
                if i == strong_user
                    % sum_rate = sum_rate + log2(1 + zeta(k)) - zeta(k) ...
                    %     + 2 * s(k) * sqrt((1 + zeta(k)) * H_c(k,i) * eta(k)) ...
                    %     - s(k)^2 * (I_tot(k,i) + H_c(k,i) * (inter_b + eta(k)));

                    (1 - eta(k)) * rho * h_c_eh(k,i)+delta_G  >= eh;
                    
                    % QoS constraint for backscatter
                    H_c(k,i) * eta(k)+delta_G >= R_min_back * (I_tot(k,i) + H_c(k,i) * inter_b);
                end
            end
        end
        
        % Power allocation constraints (sum of p per cluster = 1)
        for k = 1:K
            sum(p(k, :)) == 1;
        end
        
        % Reflection coefficient bounds
        eta+delta_G >= 0;
        eta <= 1+delta_G;
        
        % Maximize objective
        minimize(delta_G)
        
    cvx_end
    
    % Output
    obj = cvx_optval;
    p_alloc = p;  % Output renamed to p_alloc
    
    % Check feasibility
    if ~strcmp(cvx_status, 'Solved')
        fprintf('Warning: CVX status = %s\n', cvx_status);
        % Fall back to initial values if optimization failed
        if exist('p_init', 'var')
            p_alloc = p_init;
        end
        if exist('eta_init', 'var')
            eta = eta_init;
        end
    end
end