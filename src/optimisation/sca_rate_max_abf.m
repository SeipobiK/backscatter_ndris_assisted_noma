function [W_opt, A_opt, B_opt, A_c_opt, B_c_opt, obj_prev, status] = sca_rate_max_abf(...
    para,channel_data, H, H_c, A_prev, B_prev, A_c_prev, B_c_prev, decoding_order, alpha)

    % Parameters
    K   = para.K;
    K_c = para.K_c;
    M   = para.M;
    eta = para.eta;
    P_max = para.P_max;
    noise = para.noise;
    R_min = para.R_min_n;
    R_min_c = para.R_min_c;
    FT= para.FT;
    eh=para.bst_threshold;
    rho=para.rho;
    
    intra_i = zeros(K, K_c); % Initialize A_n vector
    inteer_i = zeros(K, K_c); % Initialize B_n vector
    inteer_b = zeros(K, K_c); % Initialize B_n vector
    inteer_b_all = zeros(K, K_c); % Initialize B_n vector

    cvx_begin  quiet
        cvx_solver mosek
        % Variables
        
        variable W(M,M,K) hermitian semidefinite
        variable A(K,K_c) nonnegative
        variable B(K,K_c) nonnegative
        variable A_c(K) nonnegative      % Backscatter signal variable
        variable B_c(K) nonnegative      % Backscatter interference variable
        variable R_c(K) nonnegative      % Backscatter rate variable
        % variable delta_g nonnegative
        variable R(K,K_c) nonnegative
        
        sum_rate = 0;
        for k = 1:K
            for i = 1:K_c
                sum_rate = sum_rate + R(k,i);
            end
            sum_rate = sum_rate + R_c(k); % only once
        end

        
        % Objective
        maximize(sum_rate)

        subject to

        %% ===== POWER CONSTRAINT =====
        sum_power = 0;
        for k = 1:K
            sum_power = sum_power + real(trace(W(:,:,k)));
        end
        sum_power <= P_max ;

        %% ===== MAIN LOOP =====
        for k = 1:K
            % for j = i+1:K_c
            %     R(k,i) - R(k,j) <= FT ;
            %     R(k,j) - R(k,i) <= FT ;
            % end  

            order_k = decoding_order(k,:); % weak → strong (ascending order)
            
            % Identify strong user (last user in decoding order)
            strong_user = order_k(end);  % Strong user has highest index in decoding order
            
            for i = 1:K_c

                %% ---------- SCA RATE (NOMA) ----------
                R(k,i)  <= ...
                    log2(1 + 1/(A_prev(k,i)*B_prev(k,i))) ...
                    - (log2(exp(1))/(A_prev(k,i)*(1 + A_prev(k,i)*B_prev(k,i)))) ...
                    * (A(k,i) - A_prev(k,i)) ...
                    - (log2(exp(1))/(B_prev(k,i)*(1 + A_prev(k,i)*B_prev(k,i)))) ...
                    * (B(k,i) - B_prev(k,i));
                
                %% QoS for NOMA users
                R(k,i) >= R_min ;
                
                %% ---------- INTER-CLUSTER INTERFERENCE ----------
                inter = 0;
                inter_b = 0;

                for j = 1:K
                    if j ~= k
                        inter   = inter  + real(trace(W(:,:,j) * H{k,i}'*H{k,i}));
                        inter_b = inter_b + real(trace(W(:,:,j) * H_c{k,i}'*H_c{k,i})) * eta(j);
                    end
                end

                
                
                %% ---------- INTRA-CLUSTER INTERFERENCE (NOMA) ----------
                intra = 0;
                pos = find(order_k == i);
                h_i = H{k,i};
                H_i = h_i'*h_i;

                for idx = pos+1:K_c
                    j_user = order_k(idx);
                    intra = intra + alpha(k,j_user) * real(trace(W(:,:,k) * H_i));
                    % abs(R(k,i)-R(k,j_user))<=1;
                end

                

                %% ---------- NOMA SIGNAL ----------
                signal = real(trace(W(:,:,k) * H{k,i}'*H{k,i})) * alpha(k,i);
                inv_pos(A(k,i)) <= signal ;

                %% ---------- NOMA INTERFERENCE ----------
                B(k,i) >= intra + inter + inter_b ...
                          + real(trace(W(:,:,k) * H_c{k,i}'*H_c{k,i})) * eta(k) ...
                          + noise;

                %% ===== BACKSCATTER CONSTRAINTS (ONLY FOR STRONG USER) =====
                if i == strong_user
                    eh<=(1-eta(k))*rho*real(trace(W(:,:,k) * (H_c{k,i}/channel_data.f{k,i})'*(H_c{k,i}/channel_data.f{k,i})));
                    % Backscatter signal constraint
                    signal_c = eta(k) * real(trace(W(:,:,k) * H_c{k,i}' * H_c{k,i}));
                    % disp(['user doing backscatter: ',num2str(strong_user )]);
        
                    
                    % Use SCA for backscatter rate if needed
                    R_c(k)   <= ...
                        log2(1 + 1/(A_c_prev(k)*B_c_prev(k))) ...
                        - (log2(exp(1))/(A_c_prev(k)*(1 + A_c_prev(k)*B_c_prev(k)))) ...
                        * (A_c(k) - A_c_prev(k)) ...
                        - (log2(exp(1))/(B_c_prev(k)*(1 + A_c_prev(k)*B_c_prev(k)))) ...
                        * (B_c(k) - B_c_prev(k));
                        
                    inv_pos(A_c(k)) <= signal_c ; 
                        % Backscatter interference constraint
                    B_c(k)  >= inter + inter_b + noise;
                    
                    % Backscatter QoS (only for strong user)
                    R_c(k) >= R_min_c ;
                end
                
            end
        end

    cvx_end

    %% ===== OUTPUT =====
    obj_prev = cvx_optval;
    A_opt = A;
    B_opt = B;
    A_c_opt = A_c;
    B_c_opt = B_c;
    W_opt = W;
    status = cvx_status;

    % for k=1:K
        
    %         % disp(['Near user:CVX ',num2str(R(k,1))]);
    %         % disp(['Far user:CVBX ',num2str(R(k,2))]);
    %         % disp(['Backscatter user:CVX ',num2str(R_c(k))]);
    %     % disp(R_c(k));
    % end

end
