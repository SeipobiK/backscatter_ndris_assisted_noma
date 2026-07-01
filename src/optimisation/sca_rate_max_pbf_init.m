% function [V_opt, A_opt, B_opt, A_c_opt, B_c_opt, obj_prev, status] =sca_rate_max_pbf_init(para,w_k,channel_data,decoding_order,...
%     A_prev, B_prev, A_c_prev, B_c_prev,alpha,J_r,J_t,eta)

%     %% Parameters
%     K       = para.K;
%     K_c     = para.K_c;
%     N       = para.N;
%     noise   = para.noise;
%     R_min   = para.R_min_n;
%     R_min_c = para.R_min_c;
%     eh      = para.bst_threshold;
%     rho     = para.rho;

%     %% Precompute cascaded channels
%     H   = cell(K, K_c);
%     H_c = cell(K, K_c);

%     for k = 1:K
%         for i = 1:K_c
%             H{k,i} = diag(channel_data.g{k,i}' * J_r) * ...
%                      J_t * channel_data.H_all * w_k(:,k);

%             H_c{k,i} = diag(channel_data.g_b{k}' * J_r) * ...
%                        J_t * channel_data.H_all * ...
%                        channel_data.f{k,i} * w_k(:,k);
%         end
%     end

%     %% CVX optimization
%     cvx_begin quiet
%         cvx_solver mosek

%         %% Variables
%         variable V(N,N) hermitian semidefinite
%         variable A(K,K_c) nonnegative
%         variable B(K,K_c) nonnegative
%         variable A_c(K) nonnegative
%         variable B_c(K) nonnegative
%         variable R_c(K) nonnegative
%         variable R(K,K_c) nonnegative

%         %% Dual variables
%         dual variable lambda_rank
%         dual variables lambda_diag{N}

%         dual variables lambda_rate{K,K_c}
%         dual variables lambda_qos{K,K_c}
%         dual variables lambda_signal{K,K_c}
%         dual variables lambda_interf{K,K_c}

%         dual variables lambda_eh{K}
%         dual variables lambda_c_signal{K}
%         dual variables lambda_c_rate{K}
%         dual variables lambda_c_interf{K}
%         dual variables lambda_c_qos{K}

%         %% Objective
%         sum_rate = 0;

%         for k = 1:K
%             for i = 1:K_c
%                 sum_rate = sum_rate + R(k,i);
%             end
%             sum_rate = sum_rate + R_c(k);
%         end

%         maximize(sum_rate)

%         subject to

%             %% Rank-one relaxation / penalty constraint
%             % lambda_rank : V_max' * V * V_max >= epsln_1 * trace(V);

%             %% Unit-modulus diagonal constraints
%             for m = 1:N
%                 lambda_diag{m} : V(m,m) == 1;
%             end

%             %% Main constraints
%             for k = 1:K

%                 order_k = decoding_order(k,:);   % weak -> strong
%                 strong_user = order_k(end);      % SIC/backscatter decoding user

%                 for i = 1:K_c

%                     %% ---------- SCA RATE CONSTRAINT FOR NOMA USER ----------
%                     lambda_rate{k,i} : R(k,i) <= ...
%                         log2(1 + 1/(A_prev(k,i) * B_prev(k,i))) ...
%                         - (log2(exp(1)) / ...
%                         (A_prev(k,i) * (1 + A_prev(k,i) * B_prev(k,i)))) ...
%                         * (A(k,i) - A_prev(k,i)) ...
%                         - (log2(exp(1)) / ...
%                         (B_prev(k,i) * (1 + A_prev(k,i) * B_prev(k,i)))) ...
%                         * (B(k,i) - B_prev(k,i));

%                     %% ---------- NOMA QoS CONSTRAINT ----------
%                     lambda_qos{k,i} : R(k,i) >= R_min;

%                     %% ---------- INTER-CLUSTER INTERFERENCE ----------
%                     inter = 0;
%                     inter_b = 0;

%                     for j = 1:K
%                         if j ~= k

%                             h_inter = diag(channel_data.g{k,i}' * J_r) * ...
%                                       J_t * channel_data.H_all * w_k(:,j);

%                             inter = inter + real(trace(V * (h_inter * h_inter')));

%                             h_inter_b = diag(channel_data.g_b{k}' * J_r) * ...
%                                         J_t * channel_data.H_all * ...
%                                         channel_data.f{k,i} * w_k(:,j);

%                             inter_b = inter_b + ...
%                                       real(trace(V * (h_inter_b * h_inter_b'))) ...
%                                       * eta(j);
%                         end
%                     end

%                     %% ---------- INTRA-CLUSTER NOMA INTERFERENCE ----------
%                     intra = 0;

%                     pos = find(order_k == i);

%                     h_i = H{k,i};
%                     H_i = h_i * h_i';

%                     for idx = pos+1:K_c
%                         j_user = order_k(idx);
%                         intra = intra + alpha(k,j_user) * real(trace(V * H_i));
%                     end

%                     %% ---------- NOMA SIGNAL CONSTRAINT ----------
%                     signal = alpha(k,i) * real(trace(V * H_i));

%                     lambda_signal{k,i} : inv_pos(A(k,i)) <= signal;

%                     %% ---------- NOMA INTERFERENCE CONSTRAINT ----------
%                     lambda_interf{k,i} : B(k,i) >= ...
%                         intra + inter + inter_b ...
%                         + eta(k) * real(trace(V * H_c{k,i} * H_c{k,i}')) ...
%                         + noise;

%                     %% ---------- BACKSCATTER CONSTRAINTS FOR STRONG USER ----------
%                     if i == strong_user

%                         h_eh = diag(channel_data.g_b{k}' * J_r) * ...
%                                J_t * channel_data.H_all * w_k(:,k);

%                         %% Energy harvesting constraint
%                         lambda_eh{k} : eh <= ...
%                             (1 - eta(k)) * rho * real(trace(V * (h_eh * h_eh')));

%                         %% Backscatter signal constraint
%                         signal_c = eta(k) * real(trace(V * H_c{k,i} * H_c{k,i}'));

%                         lambda_c_signal{k} : inv_pos(A_c(k)) <= signal_c;

%                         %% SCA rate constraint for backscatter
%                         lambda_c_rate{k} : R_c(k) <= ...
%                             log2(1 + 1/(A_c_prev(k) * B_c_prev(k))) ...
%                             - (log2(exp(1)) / ...
%                             (A_c_prev(k) * (1 + A_c_prev(k) * B_c_prev(k)))) ...
%                             * (A_c(k) - A_c_prev(k)) ...
%                             - (log2(exp(1)) / ...
%                             (B_c_prev(k) * (1 + A_c_prev(k) * B_c_prev(k)))) ...
%                             * (B_c(k) - B_c_prev(k));

%                         %% Backscatter interference constraint
%                         lambda_c_interf{k} : B_c(k) >= inter + inter_b + noise;

%                         %% Backscatter QoS constraint
%                         lambda_c_qos{k} : R_c(k) >= R_min_c;
%                     end
%                 end
%             end

%     cvx_end

%     %% Outputs
%     obj_prev = cvx_optval;
%     status   = cvx_status;

%     V_opt   = V;
%     A_opt   = A;
%     B_opt   = B;
%     A_c_opt = A_c;
%     B_c_opt = B_c;

%     %% Store dual variables
%     % dual_vars.lambda_rank     = lambda_rank;
%     dual_vars.lambda_diag     = lambda_diag;

%     dual_vars.lambda_rate     = lambda_rate;
%     dual_vars.lambda_qos      = lambda_qos;
%     dual_vars.lambda_signal   = lambda_signal;
%     dual_vars.lambda_interf   = lambda_interf;

%     dual_vars.lambda_eh       = lambda_eh;
%     dual_vars.lambda_c_signal = lambda_c_signal;
%     dual_vars.lambda_c_rate   = lambda_c_rate;
%     dual_vars.lambda_c_interf = lambda_c_interf;
%     dual_vars.lambda_c_qos    = lambda_c_qos;
%     % % Display dual variables for debugging  
%     % disp('Dual Variables:');
%     % % disp('lambda_rank:');
%     % % disp(lambda_rank);
%     % disp('lambda_diag:');
%     % disp(lambda_diag);
%     % disp('lambda_rate:');
%     % disp(lambda_rate);
%     % disp('lambda_qos:');
%     % disp(lambda_qos);
%     % disp('lambda_signal:');
%     % disp(lambda_signal);
%     % disp('lambda_interf:');
%     % disp(lambda_interf);
%     % disp('lambda_eh:');
%     % disp(lambda_eh);
%     % disp('lambda_c_signal:');
%     % disp(lambda_c_signal);
%     % disp('lambda_c_rate:');
%     % disp(lambda_c_rate);
%     % disp('lambda_c_interf:');
%     % disp(lambda_c_interf);
%     % disp('lambda_c_qos:');
%     % disp(lambda_c_qos); 

% end


function [V_opt, A_opt, B_opt, A_c_opt, B_c_opt, obj_prev, status] =sca_rate_max_pbf_init(para,w_k,channel_data,decoding_order,...
    A_prev, B_prev, A_c_prev, B_c_prev,alpha,J_r,J_t,eta)
    % Parameters
    K   = para.K;
    K_c = para.K_c;
    N   = para.N;
    % eta = para.eta;
    noise = para.noise;
    R_min = para.R_min_n;
    R_min_c = para.R_min_c;
    FT= para.FT;
    eh=para.bst_threshold;
    rho=para.rho;


   H= cell(K, K_c);
   H_c = cell(K, K_c);
   for k=1:K
        for i=1:K_c
            H{k,i}  = diag(channel_data.g{k,i}'*J_r)*J_t*channel_data.H_all*w_k(:, k);
            H_c{k,i}  = diag(channel_data.g_b{k}'*J_r)*J_t*channel_data.H_all * channel_data.f{k,i}*w_k(:,k);
        end
   end

   cvx_begin quiet 
        cvx_solver mosek
        % Variables
        
        variable V(N,N) hermitian semidefinite
        variable A(K,K_c) nonnegative
        variable B(K,K_c) nonnegative
        variable A_c(K) nonnegative      % Backscatter signal variable
        variable B_c(K) nonnegative      % Backscatter interference variable
        variable R_c(K) nonnegative      % Backscatter rate variable
        variable R(K,K_c) nonnegative


        % Objective
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

                for m=1:N
                    V(m,m) == 1;
                end

                for k = 1:K
                    V * eye(N) == hermitian_semidefinite(N);
                end

                %% ===== MAIN LOOP =====
                for k = 1:K

                    order_k = decoding_order(k,:); % weak → strong (ascending order)
                    
                    % Identify strong user (last user in decoding order)
                    strong_user = order_k(end);  % Strong user has highest index in decoding order
                    
                    for i = 1:K_c
                    % for j = i+1:K_c
                    %     R(k,i) - R(k,j) <= FT ;
                    %     R(k,j) - R(k,i) <= FT ;
                    % end  
                        %% ---------- SCA RATE (NOMA) ----------
                        % Constarint 3
                        R(k,i)  <= ...
                            log2(1 + 1/(A_prev(k,i)*B_prev(k,i))) ...
                            - (log2(exp(1))/(A_prev(k,i)*(1 + A_prev(k,i)*B_prev(k,i)))) ...
                            * (A(k,i) - A_prev(k,i)) ...
                            - (log2(exp(1))/(B_prev(k,i)*(1 + A_prev(k,i)*B_prev(k,i)))) ...
                            * (B(k,i) - B_prev(k,i));
                        
                        %% QoS for NOMA users
                        % Constarint 4
                        R(k,i) >= R_min;
                        
                        %% ---------- INTER-CLUSTER INTERFERENCE ----------
                        inter = 0;
                        inter_b = 0;

                        for j = 1:K
                            if j ~= k
                                inter   = inter  + real(trace(V * (diag(channel_data.g{k,i}'*J_r)*J_t*channel_data.H_all*w_k(:, j)) ...
                                *(diag(channel_data.g{k,i}'*J_r)*J_t*channel_data.H_all*w_k(:, j))'));
                                inter_b = inter_b + real(trace(V * ( diag(channel_data.g_b{k}'*J_r)*J_t*channel_data.H_all * channel_data.f{k,i}*w_k(:,j))...
                                *(diag(channel_data.g_b{k}'*J_r)*J_t*channel_data.H_all * channel_data.f{k,i}*w_k(:,j))')) * eta(k);
                            end
                        end
                        
                        

                        %% ---------- INTRA-CLUSTER INTERFERENCE (NOMA) ----------
                        intra = 0;
                        pos = find(order_k == i);
                        h_i = H{k,i};
                        H_i = h_i*h_i';

                        for idx = pos+1:K_c
                            j_user = order_k(idx);
                            intra = intra + alpha(k,j_user) * real(trace(V * H_i));
                            % abs(R(k,i)-R(k,j_user))<=1;

                            % disp(['alpha strong user  ',num2str(alpha(k,j_user) )]);
                            % disp(['index j_user user  ',num2str( j_user )]);
                            % disp(['alpha pos user  ',num2str( pos )]);
                        end
                        
                        
                        %% ---------- NOMA SIGNAL ----------

                        % Constarint 5
                        signal = real(trace(V * H_i)) * alpha(k,i);
                        inv_pos(A(k,i)) <= signal;

                        %% ---------- NOMA INTERFERENCE ----------
                        % Constarint 6
                        B(k,i) >= intra + inter + inter_b ...
                                + real(trace(V * H_c{k,i}*H_c{k,i}')) * eta(k) ...
                                + noise;

                        %% ===== BACKSCATTER CONSTRAINTS (ONLY FOR STRONG USER) =====
                        if i == strong_user
                            eh<=(1-eta(k))*rho*real(trace(V * ( diag(channel_data.g_b{k}'*J_r)*J_t*channel_data.H_all*w_k(:,k))...
                                *(diag(channel_data.g_b{k}'*J_r)*J_t*channel_data.H_all*w_k(:,k))'));
                            % Backscatter signal constraint
                            signal_c = eta(k) * real(trace(V * H_c{k,i} * H_c{k,i}'));
                            % inv_pos(A_c(k)) <= signal_c ;

                            % disp(['intra user  ',num2str( intra )]);
                            
                            % Use SCA for backscatter rate if needed
                            R_c(k)   <= ...
                                log2(1 + 1/(A_c_prev(k)*B_c_prev(k))) ...
                                - (log2(exp(1))/(A_c_prev(k)*(1 + A_c_prev(k)*B_c_prev(k)))) ...
                                * (A_c(k) - A_c_prev(k)) ...
                                - (log2(exp(1))/(B_c_prev(k)*(1 + A_c_prev(k)*B_c_prev(k)))) ...
                                * (B_c(k) - B_c_prev(k));
                                
                            inv_pos(A_c(k)) <= signal_c; 
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
    V_opt = V;
    status = cvx_status;
end