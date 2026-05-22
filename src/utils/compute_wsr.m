function [sum_rate,R,R_c,A,A_c,B,B_c,intra_i,inteer_i,inteer_b,inteer_b_all] = compute_wsr(...
    para, H, H_c,decoding_order, alpha,w_k,eta)
   % Parameters
    K   = para.K;
    K_c = para.K_c;
    M   = para.M;
    % eta = para.eta;
    P_max = para.P_max;
    noise = para.noise;

    intra_i = zeros(K, K_c); % Initialize A_n vector
    inteer_i = zeros(K, K_c); % Initialize B_n vector
    inteer_b = zeros(K, K_c); % Initialize B_n vector
    inteer_b_all = zeros(K, K_c); % Initialize B_n vector

    A = zeros(K, K_c); % Initialize A_n vector
    B = zeros(K, K_c); % Initialize B_n vector
    A_c = zeros(K,1); % Initialize A_f vector
    B_c = zeros(K,1); % Initialize B_f vector
    R=zeros(K, K_c); % Initialize R_n vector
    R_c=zeros(K,1); % Initialize R_f vector   
        
    sum_rate=0;


        %% ===== MAIN LOOP =====
        for k = 1:K

            order_k = decoding_order(k,:); % weak → strong (ascending order)
            
            % Identify strong user (last user in decoding order)
            strong_user = order_k(end);  % Strong user has highest index in decoding order
            
            for i = 1:K_c
                
                %% ---------- INTER-CLUSTER INTERFERENCE ----------
                inter = 0;
                inter_b = 0;

                for j = 1:K
                    if j ~= k
                        inter   = inter  + abs(H{k,i}*w_k(:, j)).^2;
                        inter_b = inter_b + abs(H_c{k,i}*w_k(:, j)).^2* eta(j);
                    end
                end

                %% ---------- INTRA-CLUSTER INTERFERENCE (NOMA) ----------
                intra = 0;
                pos = find(order_k == i);
                H_i = abs(H{k,i}*w_k(:, k)).^2;

                for idx = pos+1:K_c
                    j_user = order_k(idx);
                    intra = intra + alpha(k,j_user) * H_i;
                end

                 % Minimum SINR for near user

                %% ---------- NOMA SIGNAL ----------
                signal =  abs(H{k,i}*w_k(:, k)).^2 * alpha(k,i);
                A(k,i)= signal ;

                %% ---------- NOMA INTERFERENCE ----------
                B(k,i) = intra + inter + inter_b ...
                          + abs(H_c{k,i}*w_k(:, k)).^2* eta(k) ...
                          + noise;
                R(k,i) = log2(1 + A(k,i)/B(k,i));
                sum_rate = sum_rate + R(k,i);
                
                gamma_noma=2^para.R_min_n - 1;

                weak_user   = order_k(1);
                strong_user = order_k(end);

                H_weak   = abs(H{k,weak_user}   * w_k(:,k))^2;
                H_strong = abs(H{k,strong_user} * w_k(:,k))^2;

                Hc_weak   = abs(H_c{k,weak_user}   * w_k(:,k))^2;
                Hc_strong = abs(H_c{k,strong_user} * w_k(:,k))^2;

                % If inter_b already includes own BD, subtract own term to avoid double counting
                inter_b_weak_others   = inter_b - Hc_weak   * eta(k);
                inter_b_strong_others = inter_b - Hc_strong * eta(k);

                alpha_min_weak = gamma_noma * ( ...
                    H_weak + Hc_weak*eta(k) + inter_b_weak_others + inter + noise ) ...
                    / (H_weak * (1 + gamma_noma));

                alpha_min_strong = gamma_noma * ( ...
                    Hc_strong*eta(k) + inter_b_strong_others + inter + noise ) ...
                    / H_strong;

                alpha_max_weak = 1 - alpha_min_strong;
                % Since alpha_weak + alpha_strong = 1
                alpha_max_weak   = 1 - alpha_min_strong;
                alpha_max_strong = 1 - alpha_min_weak;

                % fprintf('\nCluster %d\n', k);
                % fprintf('Weak user   = %d\n', weak_user);
                % fprintf('Strong user = %d\n', strong_user);

                % fprintf('alpha_weak   range: [%.6f, %.6f]\n', ...
                %     alpha_min_weak, alpha_max_weak);

                % fprintf('alpha_strong range: [%.6f, %.6f]\n', ...
                %     alpha_min_strong, alpha_max_strong);

                fprintf('Current alpha weak   = %.6f\n', alpha(k,weak_user));
                fprintf('Current alpha strong = %.6f\n', alpha(k,strong_user));


                 alpha_ = power_allocation_opt(...
                 para, H, H_c,decoding_order,w_k);
                fprintf('Current alpha weak (from assign)  = %.6f\n', alpha_(k,weak_user));
                fprintf('Current alpha strong (from assign) = %.6f\n', alpha_(k,strong_user));

                % if alpha_min_weak <= alpha_max_weak && alpha_min_strong <= alpha_max_strong
                %     fprintf('Status: FEASIBLE\n');
                % else
                %     fprintf('Status: INFEASIBLE\n');
                % end              

                    order_k = decoding_order(k,:);
 

                
               

                %% ===== BACKSCATTER CONSTRAINTS (ONLY FOR STRONG USER) =====
                if i == strong_user

                    A_c(k) = abs(H_c{k,i}*w_k(:, k)).^2* eta(k); 
                        % Backscatter interference constraint
                    B_c(k)  = inter + inter_b + noise;
                    R_c(k) = log2(1 + A_c(k)/B_c(k));
                    sum_rate = sum_rate + R_c(k);
                     
                end
            

                    intra_i(k, i)=intra; % Initialize A_n vector
                    inteer_i(k, i) = inter + noise; % Initialize B_n vector
                    inteer_b(k, i) =inter_b; % Initialize B_n vector
                    inteer_b_all(k, i) = inter_b + abs(H_c{k,i}*w_k(:, k)).^2* eta(k); % Initialize B_n vector
            end
        end

end
