function alpha = assign_power_by_decoding_order(para, decoding_order)
    % assign_power_by_decoding_order - Fixed power allocation based on decoding order
    % Weaker user (lower decoding order, decodes first) gets MORE power
    % Stronger user (higher decoding order, decodes last) gets LESS power
    % Input:
    %   para - simulation parameters (contains K, K_c)
    %   decoding_order - decoding order for each cluster (K x K_c)
    %                     Position 1 = weakest (decodes first)
    %                     Position K_c = strongest (decodes last)
    % Output:
    %   alpha - power allocation coefficients (K x K_c)
    
    K = para.K;
    K_c = para.K_c;
    alpha = zeros(K, K_c);
    
    for k = 1:K
        order = decoding_order(k, :);
        
        if K_c == 2
            % For 2 users per cluster
            % Weakest (position 1, decodes first) gets 70%
            % Strongest (position 2, decodes last) gets 30%
            alpha(k, order(1)) = 0.5;   % Weaker user
            alpha(k, order(2)) = 0.5;   % Stronger user

            % alpha(k, 1) = 0.5;   % Weaker user
            % alpha(k, 2) = 0.5;   % Stronger user
            
        elseif K_c == 3 
            % For 3 users per cluster
            % Weakest (position 1, decodes first) gets 50%
            % Middle (position 2) gets 30%
            % Strongest (position 3, decodes last) gets 20%
            alpha(k, order(1)) = 0.35;   % Weakest user
            alpha(k, order(2)) = 0.35;   % Middle user
            alpha(k, order(3)) = 0.3;   % Strongest user

            % alpha(k, 1) = 0.35;   % Weakest user
            % alpha(k, 2) = 0.35;   % Middle user
            % alpha(k, 3) = 0.3;   % Strongest user
        else
            % Default fallback for other K_c values
            % Linear decreasing: weakest gets most, strongest gets least
            total_power = 1;
            power_weights = K_c:-1:1;
            power_weights = power_weights / sum(power_weights);
            for idx = 1:K_c
                user_idx = order(idx);
                alpha(k, user_idx) = total_power * power_weights(idx);
            end

            alpha(k, 1) = 0.1;   % Weakest user
            alpha(k, 2) = 0.2;   % Middle user
            alpha(k, 3) = 0.3;   % Middle user
            alpha(k, 4) = 0.4;   % Strongest user
        end
    end
    
    % Ensure sum of power = 1 for each cluster (numerical safety)
    for k = 1:K
        if abs(sum(alpha(k, :)) - 1) > 1e-6
            alpha(k, :) = alpha(k, :) / sum(alpha(k, :));
        end
    end
end