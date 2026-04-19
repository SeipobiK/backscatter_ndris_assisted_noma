function alpha = power_allocation_opt(para, H, H_c, decoding_order, w_k)
    % Power allocation based on decoding order (SIC)
    % Weaker users (lower decoding order) get minimum power to meet QoS
    % Stronger users (higher decoding order) get remaining power
    
    % Parameters
    K = para.K;
    K_c = para.K_c;
    eta = para.eta;
    noise = para.noise;
    % R_min = para.R_min_n;  % Minimum rate requirement
    R_min=2^para.R_min_n-1;
    
    alpha = zeros(K, K_c);
    effective_gain = zeros(K, K_c);
    
    for k = 1:K
        order = decoding_order(k, :);  % order(1) = weakest, order(end) = strongest
        
        %% Calculate effective channel gain for each user
        for idx = 1:K_c
            user = order(idx);
            
            % Inter-cluster interference
            inter = 0;
            inter_b = 0;
            
            for j = 1:K
                if j ~= k
                    inter = inter + abs(H{k, user} * w_k(:, j))^2;
                    inter_b = inter_b + abs(H_c{k, user} * w_k(:, j))^2 * eta(k);
                end
            end
            
            % Backscatter interference from own cluster
            backscatter = abs(H_c{k, user} * w_k(:, k))^2 * eta(k);
            
            % Total interference + noise
            total_interference = inter + inter_b + backscatter ;


            
            % Effective channel gain: |h|^2 / (I + noise)
            effective_gain(k, user) = abs(H{k, user} * w_k(:, k))^2 / total_interference;
        end
        
        % Weakest user (decodes first) - gets minimum power to meet QoS
        weakest_user = order(1);
        alpha(k, weakest_user) = (R_min / (1 + R_min)) * (1 + 1/effective_gain(k, weakest_user));
        
        % % Clamp to reasonable range
        % alpha(k, weakest_user) = max(alpha(k, weakest_user), 0.05);
        % alpha(k, weakest_user) = min(alpha(k, weakest_user), 0.5);
        
        % Remaining users get remaining power
        remaining_power = 1 - alpha(k, weakest_user);
        
        if K_c == 2
            % 2 users: strongest gets all remaining
            strongest_user = order(2);
            alpha(k, strongest_user) = remaining_power;
            
        elseif K_c == 3
            % 3 users: allocate based on gains
            middle_user = order(2);
            strongest_user = order(3);
            
            % Middle user gets portion based on its gain
            gain_mid = effective_gain(k, middle_user);
            gain_strong = effective_gain(k, strongest_user);
            total_gain = gain_mid + gain_strong;
            
            alpha(k, middle_user) = remaining_power * (gain_mid / total_gain);
            alpha(k, strongest_user) = remaining_power * (gain_strong / total_gain);
            
        else
            % General case: allocate proportionally to gains
            total_gain = 0;
            for idx = 2:K_c
                total_gain = total_gain + effective_gain(k, order(idx));
            end
            
            for idx = 2:K_c
                user = order(idx);
                alpha(k, user) = remaining_power * (effective_gain(k, user) / total_gain);
            end
        end
        
        % Normalize to ensure sum = 1
        alpha(k, :) = alpha(k, :) / sum(alpha(k, :));
        
        % Verify QoS for weakest user
        actual_sinr = effective_gain(k, weakest_user) * alpha(k, weakest_user);
        target_sinr = 2^R_min - 1;
        
        if actual_sinr < target_sinr - 1e-6
            fprintf('Warning: Cluster %d, User %d SINR = %.4f < target %.4f\n', ...
                k, weakest_user, actual_sinr, target_sinr);
        end
    end
end