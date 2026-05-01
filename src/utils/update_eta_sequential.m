function eta_new = update_eta_sequential(para, eta_old, y, z, zeta, h_ki_c, h_ji_c, A_ki, I_tot, I_btot, i_star)
% UPDATE_ETA_SEQUENTIAL - Sequential update for reflection coefficients
%
% Computes gamma_min from rate requirements:
%   gamma_min_ki = 2^para.R_min_n - 1  (for NOMA users)
%   gamma_min_c  = 2^para.R_min_c - 1  (for backscatter)

    K = para.K;
    K_c = para.K_c;
    
    % ===== COMPUTE SINR REQUIREMENTS FROM RATE REQUIREMENTS =====
    % For NOMA users (cellular links)
    if isfield(para, 'R_min_n')
        gamma_min_ki = 2^para.R_min_n - 1;
        % If R_min_n is a scalar, expand to matrix
        if isscalar(gamma_min_ki)
            gamma_min_ki = gamma_min_ki * ones(K, K_c);
        end
    else
        % Default if not provided
        gamma_min_ki = 0.1 * ones(K, K_c);
    end
    
    % For backscatter link
    if isfield(para, 'R_min_c')
        gamma_min_c = 2^para.R_min_c - 1;
        if isscalar(gamma_min_c)
            gamma_min_c = gamma_min_c * ones(K, 1);
        end
    else
        % Default if not provided
        gamma_min_c = 0.1 * ones(K, 1);
    end
    
    eta_new = eta_old;
    
    % ===== UPDATE EACH USER SEQUENTIALLY =====
    for k = 1:K
        
        i_st = i_star(k);
        
        % -------------------------------------------------------------
        % STEP 1: Compute C_k = sum_i y_{k,i}^2 * h_{k,i}^c + z_k^2 * h_{k,i*}^c
        % -------------------------------------------------------------
        sum_y2_h = 0;
        for i = 1:K_c
            sum_y2_h = sum_y2_h + y(k, i)^2 * h_ki_c(k, i);
        end
        
        C_k = sum_y2_h + z(k)^2 * h_ki_c(k, i_st);
        
        if C_k < 1e-12
            C_k = 1e-12;
        end
        
        % -------------------------------------------------------------
        % STEP 2: Unconstrained optimum (Eq. 16)
        % eta_unc = z_k^2 * (1+zeta_k) * h_{k,i*}^c / C_k^2
        % -------------------------------------------------------------
        numerator_unc = z(k)^2 * (1 + zeta(k)) * h_ki_c(k, i_st);
        eta_unc = numerator_unc / (C_k * C_k);
        
        % -------------------------------------------------------------
        % STEP 3: Lower bound from constraint C2 (backscatter QoS)
        % eta_min = gamma_min_c * (sum_{j≠k} h_{j,i*}^c * eta_j + I_btot) / h_{k,i*}^c
        % -------------------------------------------------------------
        sum_interf = 0;
        for j = 1:K
            if j ~= k
                if j < k
                    sum_interf = sum_interf + h_ji_c(j, i_st) * eta_new(j);
                else
                    sum_interf = sum_interf + h_ji_c(j, i_st) * eta_old(j);
                end
            end
        end
        
        if h_ki_c(k, i_st) > 1e-12
            eta_min = gamma_min_c(k) * (sum_interf + I_btot(k)) / h_ki_c(k, i_st);
        else
            eta_min = 0;
        end
        eta_min = max(0, eta_min);
        
        % -------------------------------------------------------------
        % STEP 4: Upper bound from constraint C1 (NOMA users' QoS)
        % For each subcarrier i:
        % eta_max = min_i [ (h_{k,i}*alpha - gamma_min_ki * I_tot) / (gamma_min_ki * h_{k,i}^c)
        %                     - sum_{j≠k} (h_{j,i}^c / h_{k,i}^c) * eta_j ]
        % -------------------------------------------------------------
        eta_max = 1;  % start with physical upper bound
        
        for i = 1:K_c
            
            if h_ki_c(k, i) > 1e-12 && gamma_min_ki(k, i) > 1e-12
                
                % Compute sum_{j≠k} (h_{j,i}^c / h_{k,i}^c) * eta_j
                sum_other = 0;
                for j = 1:K
                    if j ~= k
                        if j < k
                            sum_other = sum_other + (h_ji_c(j, i) / h_ki_c(k, i)) * eta_new(j);
                        else
                            sum_other = sum_other + (h_ji_c(j, i) / h_ki_c(k, i)) * eta_old(j);
                        end
                    end
                end
                
                % Bound from this subcarrier
                numerator_bound = A_ki(k, i) - gamma_min_ki(k, i) * I_tot(k, i);
                denominator_bound = gamma_min_ki(k, i) * h_ki_c(k, i);
                
                if denominator_bound > 1e-12
                    bound_i = numerator_bound / denominator_bound - sum_other;
                    
                    % Take the minimum bound across all subcarriers
                    if bound_i < eta_max
                        eta_max = bound_i;
                    end
                end
            end
        end
        
        % Clamp to physical bounds [0, 1]
        eta_max = min(1, max(0, eta_max));
        
        % -------------------------------------------------------------
        % STEP 5: Handle infeasible bounds (if min > max)
        % -------------------------------------------------------------
        if eta_min > eta_max + 1e-10
            % Take midpoint with small relaxation
            eta_mid = (eta_min + eta_max) / 2;
            eta_min = max(0, eta_mid - 0.1);
            eta_max = min(1, eta_mid + 0.1);
            if eta_min > eta_max
                eta_min = 0;
                eta_max = 1;
            end
        end
        
        % -------------------------------------------------------------
        % STEP 6: Projected update (Eq. 18)
        % eta_new = min(eta_max, max(eta_min, eta_unc))
        % -------------------------------------------------------------
        eta_new(k) = eta_unc;
        
        if eta_new(k) < eta_min
            eta_new(k) = eta_min;
        end
        
        if eta_new(k) > eta_max
            eta_new(k) = eta_max;
        end
        
        % Final safety clamp
        eta_new(k) = max(0, min(1, eta_new(k)));
        
    end  % end for k
    
end