function [sum_rate,R,R_c,A,A_c,B,B_c,intra_i,inteer_i,inteer_b,inteer_b_all] = compute_wsr(...
    para, H, H_c,decoding_order, alpha,w_k)
   % Parameters
    K   = para.K;
    K_c = para.K_c;
    M   = para.M;
    eta = para.eta;
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

                %% ---------- NOMA SIGNAL ----------
                signal =  abs(H{k,i}*w_k(:, k)).^2 * alpha(k,i);
                A(k,i)= signal ;

                %% ---------- NOMA INTERFERENCE ----------
                B(k,i) = intra + inter + inter_b ...
                          + abs(H_c{k,i}*w_k(:, k)).^2* eta(k) ...
                          + noise;
                R(k,i) = log2(1 + A(k,i)/B(k,i));
                sum_rate = sum_rate + R(k,i);

                %% ===== BACKSCATTER CONSTRAINTS (ONLY FOR STRONG USER) =====
                if i == strong_user

                    A_c(k) = abs(H_c{k,i}*w_k(:, k)).^2* eta(k); 
                        % Backscatter interference constraint
                    B_c(k)  = inter + inter_b + noise;
                    R_c(k) = log2(1 + A_c(k)/B_c(k));
                    sum_rate = sum_rate + R_c(k);
                end

                    intra_i(k, i)=intra; % Initialize A_n vector
                    inteer_i(k, i) = inter + noise + inter_b + abs(H_c{k,i}*w_k(:, k)).^2* eta(k); % Initialize B_n vector
                    inteer_b(k, i) = A(k,i)/B(k,i); % Initialize B_n vector
                    inteer_b_all(k, i) = inter_b + abs(H_c{k,i}*w_k(:, k)).^2* eta(k); % Initialize B_n vector
            end
        end

end
