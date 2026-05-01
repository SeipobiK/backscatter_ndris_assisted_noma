function [beta, zeta, y, z, I_tot, I_btot, h_ki, h_ki_c, h_ji_c, A_ki, i_star] = update_auxiliary_variables(...
    para, H, H_c, decoding_order, alpha, w_k, eta)

    K = para.K;
    K_c = para.K_c;
    noise = para.noise;
    
    % Initialize outputs
    beta = zeros(K, K_c);
    zeta = zeros(K, 1);
    y = zeros(K, K_c);
    z = zeros(K, 1);
    I_tot = zeros(K, K_c);
    I_btot = zeros(K, 1);
    h_ki = zeros(K, K_c);
    h_ki_c = zeros(K, K_c);
    h_ji_c = zeros(K, K_c);  % Note: this will store for current k,i but j varies
    A_ki = zeros(K, K_c);
    i_star = zeros(K, 1);
    
    % For storing interference from all users (needed for h_ji_c)
    % We'll build a 3D-like structure using cells
    h_ji_c_all = cell(K, K_c);
    
    %% ===== MAIN LOOP =====
    for k = 1:K
        order_k = decoding_order(k, :); % weak → strong (ascending order)
        
        % Identify strong user (last user in decoding order) - this is i*
        i_star(k) = order_k(end);
        
        for i = 1:K_c
            
            %% ---------- DESIRED SIGNAL POWERS ----------
            % Cellular link gain |H{k,i} * w_k|^2
            h_ki(k, i) = abs(H{k,i} * w_k(:, k))^2;
            
            % Backscatter link gain |H_c{k,i} * w_k|^2
            h_ki_c(k, i) = abs(H_c{k,i} * w_k(:, k))^2;
            
            % A_{k,i} = h_{k,i} * alpha_{k,i}
            A_ki(k, i) = h_ki(k, i) * alpha(k, i);
            
            %% ---------- INTERFERENCE COMPUTATIONS ----------
            inter_cell = 0;      % interference from other users' cellular signals
            inter_back = 0;      % interference from other users' backscatter signals
            
            for j = 1:K
                if j ~= k
                    % Cellular interference from user j
                    inter_cell = inter_cell + abs(H{k,i} * w_k(:, j))^2;
                    
                    % Backscatter interference channel gain (store for later use in h_ji_c)
                    h_ji_c_temp = abs(H_c{k,i} * w_k(:, j))^2;
                    h_ji_c_all{k, i}(j) = h_ji_c_temp;
                    
                    % Backscatter interference (with eta)
                    inter_back = inter_back + h_ji_c_temp * eta(j);
                end
            end
            
            %% ---------- INTRA-CLUSTER INTERFERENCE (NOMA) ----------
            intra = 0;
            pos = find(order_k == i);
            H_i_power = h_ki(k, i);  % reuse computed value
            
            % Users with higher decoding order (weaker users) cause intra-cluster interference
            for idx = pos+1:K_c
                j_user = order_k(idx);
                intra = intra + alpha(k, j_user) * H_i_power;
            end
            
            %% ---------- BACKSCATTER SELF-INTERFERENCE ----------
            back_self = h_ki_c(k, i) * eta(k);
            
            %% ---------- TOTAL INTERFERENCE + NOISE FOR CELLULAR LINK ----------
            I_tot(k, i) = intra + inter_cell + inter_back + back_self + noise;
            
            %% ===== LDT: UPDATE beta_{k,i} (Eq. 5) =====
            % beta_{k,i}^* = A_{k,i} / (sum_j h_{j,i}^c * eta_j + I_tot)
            % Note: sum_j h_{j,i}^c * eta_j = inter_back + back_self
            denominator_beta = inter_back + back_self + noise;
            if denominator_beta > eps
                beta(k, i) = A_ki(k, i) / denominator_beta;
            else
                beta(k, i) = A_ki(k, i) / eps;
            end
            
            %% ===== QT: UPDATE y_{k,i} (Eq. 11) =====
            % y_{k,i}^* = sqrt((1+beta_{k,i}) * A_{k,i}) / (A_{k,i} + sum_j h_{j,i}^c*eta_j + I_tot)
            denominator_y = A_ki(k, i) + inter_back + back_self + noise;
            numerator_y = sqrt(max(0, (1 + beta(k, i)) * A_ki(k, i)));
            if denominator_y > eps
                y(k, i) = numerator_y / denominator_y;
            else
                y(k, i) = 0;
            end
        end
        
        %% ===== BACKSCATTER LINK COMPUTATIONS (at i*) =====
        i_st = i_star(k);
        
        % Compute interference from other users for backscatter link at i*
        inter_back_total = 0;
        for j = 1:K
            if j ~= k
                inter_back_total = inter_back_total + abs(H_c{k, i_st} * w_k(:, j))^2 * eta(j);
            end
        end
        
        % Total interference + noise for backscatter link (denominator of SINR_c)
        I_btot(k) = inter_back_total + noise;
        
        % Backscatter signal power
        back_signal = h_ki_c(k, i_st) * eta(k);
        
        %% ===== LDT: UPDATE zeta_k (Eq. 6) =====
        % zeta_k^* = (h_{k,i*}^c * eta_k) / (sum_{j≠k} h_{j,i*}^c * eta_j + I_btot)
        denominator_zeta = inter_back_total + noise;
        if denominator_zeta > eps
            zeta(k) = back_signal / denominator_zeta;
        else
            zeta(k) = back_signal / eps;
        end
        
        %% ===== QT: UPDATE z_k (Eq. 12) =====
        % z_k^* = sqrt((1+zeta_k) * h_{k,i*}^c * eta_k) / (h_{k,i*}^c*eta_k + sum_{j≠k} h_{j,i*}^c*eta_j + I_btot)
        denominator_z = back_signal + inter_back_total + noise;
        numerator_z = sqrt(max(0, (1 + zeta(k)) * h_ki_c(k, i_st) * eta(k)));
        if denominator_z > eps
            z(k) = numerator_z / denominator_z;
        else
            z(k) = 0;
        end
    end
    
    % Reconstruct h_ji_c matrix (interference from user j to user k's backscatter link)
    % For each k,i, h_ji_c(j,i) is the interference from user j
    h_ji_c = zeros(K, K_c);
    for k = 1:K
        for i = 1:K_c
            if ~isempty(h_ji_c_all{k, i})
                % For a given (k,i), we need the interference from a specific j
                % This is used in the eta update loop. We'll store the full matrix
                for j = 1:K
                    if j ~= k && j <= length(h_ji_c_all{k, i})
                        h_ji_c(j, i) = h_ji_c_all{k, i}(j);
                    end
                end
            end
        end
    end
end