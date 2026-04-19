function [beta, t, zeta, s] = update_aux_variables_pac_rc(...
    para, H, H_c, I_tot, power, eta, decoding_order)
    
    % Parameters
    K = para.K;
    K_c = para.K_c;
    
    % Initialize
    beta = zeros(K, K_c);
    t = zeros(K, K_c);
    zeta = zeros(K, 1);
    s = zeros(K, 1);
    
    for k = 1:K
        order_k = decoding_order(k, :);
        strong_user = order_k(end);
        
        for i = 1:K_c
            % Inter-backscatter interference
            inter_b = sum(eta) - eta(k);
            
            % Intra-cluster interference from stronger users
            pos = find(order_k == i);
            intra = 0;
            for idx = pos+1:K_c
                j_user = order_k(idx);
                intra = intra + power(k, j_user);
            end
            
            % Update beta (Lagrangian dual for NOMA)
            denominator_beta = I_tot(k,i) + H_c(k,i) * (inter_b + eta(k)) + H(k,i) * intra;
            if denominator_beta > 0
                beta(k,i) = (H(k,i) * power(k,i)) / denominator_beta;
            else
                beta(k,i) = 0;
            end
            
            % Update t (quadratic transform for NOMA)
            denominator_t = I_tot(k,i) + H_c(k,i) * (inter_b + eta(k)) + H(k,i) * (intra + power(k,i));
            if denominator_t > 0
                t(k,i) = sqrt((1 + beta(k,i)) * H(k,i) * power(k,i)) / denominator_t;
            else
                t(k,i) = 0;
            end
            
            % Backscatter updates (only for strongest user)
            if i == strong_user
                % Update zeta (Lagrangian dual for backscatter)
                denominator_zeta = I_tot(k,i) + H_c(k,i) * inter_b;
                if denominator_zeta > 0
                    zeta(k) = (H_c(k,i) * eta(k)) / denominator_zeta;
                else
                    zeta(k) = 0;
                end
                
                % Update s (quadratic transform for backscatter)
                denominator_s = I_tot(k,i) + H_c(k,i) * (inter_b + eta(k));
                if denominator_s > 0
                    s(k) = sqrt((1 + zeta(k)) * H_c(k,i) * eta(k)) / denominator_s;
                else
                    s(k) = 0;
                end
            end
        end
    end
end