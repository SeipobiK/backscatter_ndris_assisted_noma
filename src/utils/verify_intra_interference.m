function verify_intra_interference(decoding_order, noma_signal, alpha, intra_i)
    % verify_intra_interference - Verify intra-cluster interference calculation
    % Input:
    %   decoding_order - decoding order (K x K_c)
    %   noma_signal - received signal power (includes alpha) = |h_eff|^2 * alpha(k,i)
    %   alpha - power allocation coefficients (K x K_c)
    %   intra_i - computed intra-cluster interference from optimization
    
    K = size(decoding_order, 1);
    K_c = size(decoding_order, 2);
    
    fprintf('\n========== INTRA-CLUSTER INTERFERENCE VERIFICATION ==========\n\n');
    
    for k = 1:K
        fprintf('Cluster %d:\n', k);
        order = decoding_order(k, :);
        
        % Calculate effective channel gain for each user
        % noma_signal = |h_eff|^2 * alpha(k,i)
        % Therefore: |h_eff|^2 = noma_signal / alpha(k,i)
        h_eff_sq = zeros(1, K_c);
        for i = 1:K_c
            if alpha(k, i) > 0
                h_eff_sq(i) = noma_signal(k, i) / alpha(k, i);
            else
                h_eff_sq(i) = 0;
            end
        end
        
        fprintf('\nEffective channel gains (|h_eff|^2):\n');
        for i = 1:K_c
            fprintf('  User %d: %.6f\n', i, h_eff_sq(i));
        end
        
        % Verify intra-cluster interference for each user
        fprintf('\nIntra-cluster interference verification:\n');
        fprintf('------------------------------------------------\n');
        
        for pos = 1:K_c
            current_user = order(pos);
            
            % Users that cause interference (higher decoding order = stronger users)
            interfering_users = order(pos+1:end);
            
            % Calculate expected intra interference
            % Intra = |h_eff|^2 * sum(alpha of interfering users)
            expected_intra = 0;
            sum_alpha_interfering = 0;
            alpha_sum_str = '';
            
            if ~isempty(interfering_users)
                alpha_values = zeros(1, length(interfering_users));
                for j = 1:length(interfering_users)
                    user_idx = interfering_users(j);
                    alpha_values(j) = alpha(k, user_idx);
                    sum_alpha_interfering = sum_alpha_interfering + alpha_values(j);
                end
                expected_intra = h_eff_sq(current_user) * sum_alpha_interfering;
                
                % Build the sum string
                alpha_sum_str = sprintf('%.3f', alpha_values(1));
                for j = 2:length(alpha_values)
                    alpha_sum_str = sprintf('%s + %.3f', alpha_sum_str, alpha_values(j));
                end
            end
            
            % Get computed intra interference
            computed_intra = intra_i(k, current_user);
            
            % Display results
            fprintf('\nUser %d (decoding position %d):\n', current_user, pos);
            fprintf('  |h_eff|^2 = %.6f\n', h_eff_sq(current_user));
            
            if ~isempty(interfering_users)
                fprintf('  Interfering users (higher decoding order): ');
                for j = 1:length(interfering_users)
                    fprintf('User %d ', interfering_users(j));
                end
                fprintf('\n');
                fprintf('  Sum of their alpha coefficients: %.6f (= %s)\n', sum_alpha_interfering, alpha_sum_str);
                fprintf('  Expected intra interference: %.6f\n', expected_intra);
            else
                fprintf('  No interfering users (weakest user)\n');
                fprintf('  Expected intra interference: 0\n');
            end
            
            fprintf('  Computed intra interference: %.6f\n', computed_intra);
            
            % Check if they match
            if abs(expected_intra - computed_intra) < 1e-6
                fprintf('  ✓ VERIFIED: Intra interference matches!\n');
            else
                fprintf('  ✗ MISMATCH: Difference = %.6f\n', abs(expected_intra - computed_intra));
            end
        end
        
        fprintf('\n========================================\n\n');
    end
end