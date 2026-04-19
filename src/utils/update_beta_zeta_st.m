function [beta,I_tot,t,h,h_c,h_eh] = update_beta_zeta_st(...
    para, H, H_c, decoding_order, alpha, w_k, eta)
    
    % Parameters
    K = para.K;
    K_c = para.K_c;
    noise = para.noise;

    
    % Initialize
    beta = zeros(K, K_c);
    t = zeros(K, K_c);
    A = zeros(K, K_c); % Initialize A_n vector
    B = zeros(K, K_c); % Initialize B_n vector

    h = zeros(K, K_c);
    h_c = zeros(K, K_c);
    h_eh=zeros(K, K_c);
    I_tot=zeros(K, K_c);
    
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
                
                beta(k,i)=A(k,i)/B(k,i);

                t(k,i)=sqrt((1+beta(k,i))*real(A(k,i)))/(A(k,i)+B(k,i));

                I_tot(k,i)=inter + noise + inter_b + abs(H_c{k,i}*w_k(:, k)).^2* eta(k);
                h(k,i)=abs(H{k,i}*w_k(:, k)).^2;
                h_c(k,i)=abs(H_c{k,i}*w_k(:, k)).^2;
                h_eh(k,i)=abs(H{k,i}*w_k(:, k)).^2;

            end
        end

end
