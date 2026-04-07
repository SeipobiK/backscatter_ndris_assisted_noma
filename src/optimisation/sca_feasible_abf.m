function [W_opt, A_opt, B_opt, A_c_opt, B_c_opt, obj_prev, status] = sca_feasible_abf(...
    para, H, H_c, A_prev, B_prev, A_c_prev, B_c_prev, decoding_order, alpha)

    % Parameters
    K   = para.K;
    K_c = para.K_c;
    M   = para.M;
    eta = para.eta;
    P_max = para.P_max;
    noise = para.noise;
    R_min = para.R_min_n;
    R_min_c = para.R_min_c;

    cvx_begin   quiet
        cvx_solver mosek
        % Variables
        
        variable W(M,M,K) hermitian semidefinite
        variable A(K,K_c) nonnegative
        variable B(K,K_c) nonnegative
        variable A_c(K) nonnegative      % Backscatter signal variable
        variable B_c(K) nonnegative      % Backscatter interference variable
        variable R_c(K) nonnegative      % Backscatter rate variable
        variable delta_g nonnegative
        variable R(K,K_c) nonnegative


        % Objective
        minimize(delta_g)

        subject to

                %% ===== POWER CONSTRAINT =====
                sum_power = 0;
                % Constarint 1
                for k = 1:K
                    sum_power = sum_power + real(trace(W(:,:,k)));
                end
                sum_power <= P_max + delta_g;
                % Constarint 2
                for k = 1:K
                    W(:, :, k) + delta_g * eye(M) == hermitian_semidefinite(M);
                end

                %% ===== MAIN LOOP =====
                for k = 1:K

                    order_k = decoding_order(k,:); % weak → strong (ascending order)
                    
                    % Identify strong user (last user in decoding order)
                    strong_user = order_k(end);  % Strong user has highest index in decoding order
                    
                    for i = 1:K_c

                        %% ---------- SCA RATE (NOMA) ----------
                        % Constarint 3
                        R(k,i) - delta_g <= ...
                            log2(1 + 1/(A_prev(k,i)*B_prev(k,i))) ...
                            - (log2(exp(1))/(A_prev(k,i)*(1 + A_prev(k,i)*B_prev(k,i)))) ...
                            * (A(k,i) - A_prev(k,i)) ...
                            - (log2(exp(1))/(B_prev(k,i)*(1 + A_prev(k,i)*B_prev(k,i)))) ...
                            * (B(k,i) - B_prev(k,i));
                        
                        %% QoS for NOMA users
                        % Constarint 4
                        R(k,i) >= R_min - delta_g;
                        
                        %% ---------- INTER-CLUSTER INTERFERENCE ----------
                        inter = 0;
                        inter_b = 0;

                        for j = 1:K
                            if j ~= k
                                inter   = inter  + real(trace(W(:,:,j) * H{k,i}'*H{k,i}));
                                inter_b = inter_b + real(trace(W(:,:,j) * H_c{k,i}'*H_c{k,i})) * eta;
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
                            % abs(R(k,i)-R(k,j_user))<=1+ delta_g;

                            % disp(['alpha strong user  ',num2str(alpha(k,j_user) )]);
                            % disp(['index j_user user  ',num2str( j_user )]);
                            % disp(['alpha pos user  ',num2str( pos )]);
                        end
                        
                        
                        %% ---------- NOMA SIGNAL ----------

                        % Constarint 5
                        signal = real(trace(W(:,:,k) * H{k,i}'*H{k,i})) * alpha(k,i);
                        inv_pos(A(k,i)) <= signal + delta_g;

                        %% ---------- NOMA INTERFERENCE ----------
                        % Constarint 6
                        B(k,i)+delta_g >= intra + inter + inter_b ...
                                + real(trace(W(:,:,k) * H_c{k,i}'*H_c{k,i})) * eta ...
                                + noise;

                        %% ===== BACKSCATTER CONSTRAINTS (ONLY FOR STRONG USER) =====
                        if i == strong_user
                            % Backscatter signal constraint
                            signal_c = eta * real(trace(W(:,:,k) * H_c{k,i}' * H_c{k,i}));
                            inv_pos(A_c(k)) <= signal_c + delta_g;

                            % disp(['intra user  ',num2str( intra )]);
                            
                            % Use SCA for backscatter rate if needed
                            R_c(k) - delta_g  <= ...
                                log2(1 + 1/(A_c_prev(k)*B_c_prev(k))) ...
                                - (log2(exp(1))/(A_c_prev(k)*(1 + A_c_prev(k)*B_c_prev(k)))) ...
                                * (A_c(k) - A_c_prev(k)) ...
                                - (log2(exp(1))/(B_c_prev(k)*(1 + A_c_prev(k)*B_c_prev(k)))) ...
                                * (B_c(k) - B_c_prev(k));
                                
                            inv_pos(A_c(k)) <= signal_c + delta_g; 
                                % Backscatter interference constraint
                            B_c(k) + delta_g >= inter + inter_b + noise;
                            
                            % Backscatter QoS (only for strong user)
                            R_c(k) >= R_min - delta_g;
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

end
