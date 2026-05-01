function [beta, zeta, y, z, sum_rate_const_alpha, sum_rate_const_eta,I_tot] = get_constants_ldt_qt1(...
    para, H, H_c, decoding_order, alpha, w_k, eta)
% =========================================================
%  Computes LDT and QT auxiliary variables:
%    beta(k,i), y(k,i)  -- for NOMA rate R_{k,i}
%    zeta(k),   z(k)    -- for BD rate   R_{k,i*}^c
% ========================================================= 
 
    %% ---- Parameters ----
    K     = para.K;
    K_c   = para.K_c;
    noise = para.noise;
 
    %% ---- Initialise ----
    A     = zeros(K, K_c);   % NOMA signal power
    B     = zeros(K, K_c);   % NOMA interference + noise
    A_c   = zeros(K, 1);     % BD signal power
    B_c   = zeros(K, 1);     % BD interference + noise
 
    beta  = zeros(K, K_c);
    I_tot  = zeros(K, K_c);
    y     = zeros(K, K_c);
    zeta  = zeros(K, 1);
    z     = zeros(K, 1);

    sum_rate_const_alpha=0;
    sum_rate_const_eta=0;
    
    %% ===== MAIN LOOP =====
    for k = 1:K
 
        order_k     = decoding_order(k,:);   % weak -> strong
        strong_user = order_k(end);          % i* for cluster k
 
        for i = 1:K_c
 
            %% ---- Inter-cluster interference ----
            % inter   : from other clusters' beamformers (no BD)
            % inter_b : BD interference from ALL BDs (j=1..K) at NOMA user
            %           = sum_{j=1}^{K} |H_c{k,i}*w_k(:,j)|^2 * eta(j)
            %           includes j=k term as well
            inter   = 0;   % from other BSs/clusters beamformers
            inter_b = 0;   % from ALL BDs -> NOMA user (k,i)
            
 
            for j = 1:K
                if j ~= k
                    % Interference from other clusters' beamformers
                    inter = inter + abs(H{k,i} * w_k(:,j)).^2;
                end
                % BD interference at NOMA user: ALL j=1..K
                inter_b = inter_b + abs(H_c{k,i} * w_k(:,j)).^2 * eta(j);
            end
 
            %% ---- Intra-cluster interference (NOMA SIC) ----
            % Users decoded after i (stronger users) cause interference
            intra = 0;
            pos   = find(order_k == i);
            H_ki  = abs(H{k,i} * w_k(:,k)).^2;   % h_{k,i}
 
            for idx = pos+1:K_c
                j_user = order_k(idx);
                intra  = intra + alpha(k, j_user) * H_ki;
            end
     
            %% ---- NOMA signal and total interference ----
            % A(k,i) = h_{k,i} * alpha_{k,i}   (signal term)
            % B(k,i) = intra + inter + inter_b  (all interference + noise)
            %        = sum_{k=1}^{K} h_{k,i}^c * eta_k  +  I_{k,i}^{tot}
            A(k,i) = H_ki * alpha(k,i);
            B(k,i) = intra + inter + inter_b + noise;

            I_tot(k,i) = inter + noise;
 
            %% ---- LDT auxiliary: beta_{k,i} = A/B (SINR) ----
            beta(k,i) = A(k,i) / B(k,i);
            sum_rate_const_alpha=sum_rate_const_alpha+log2(1+beta(k,i))-beta(k,i);
 
            %% ---- QT auxiliary: y_{k,i} ----
            % y* = sqrt((1+beta)*A) / (A + B)
            y(k,i) = sqrt((1 + beta(k,i)) * A(k,i)) / (A(k,i) + B(k,i));
 
            %% ===== BD (BACKSCATTER) TERMS — strong user i* only =====
            if i == strong_user
 
                H_c_kk = abs(H_c{k,i} * w_k(:,k)).^2;  % h_{k,i*}^c
 
                %% ---- BD signal power ----
                % A_c(k) = h_{k,i*}^c * eta_k
                A_c(k) = H_c_kk * eta(k);
 
                %% ---- BD interference from OTHER BDs only (j ~= k) ----
                % B_c(k) = sum_{j~=k} h_{j,i*}^c * eta_j + I_{k,i*}^{btot}
                % Note: inter_b already sums ALL j=1..K so subtract j=k term
                inter_b_others = inter_b - H_c_kk * eta(k);
 
                % Add non-BD interference (other beamformers + noise)
                B_c(k) = inter_b_others + inter + noise;
 
                %% ---- LDT auxiliary: zeta_k = A_c / B_c (BD SINR) ----
                zeta(k) = A_c(k) / B_c(k);
                sum_rate_const_eta=sum_rate_const_eta+log2(1+zeta(k))-zeta(k);
                %% ---- QT auxiliary: z_k ----
                % z* = sqrt((1+zeta)*A_c) / (A_c + B_c)
                z(k) = sqrt((1 + zeta(k)) * A_c(k)) / (A_c(k) + B_c(k));
            end
 
        end % i loop

                   weak = order_k(1);
            strong = order_k(end);

            x = alpha(k,weak);

            Hweak = abs(H{k,weak}*w_k(:,k))^2;
            Hstrong = abs(H{k,strong}*w_k(:,k))^2;

            a_w = 2*y(k,weak)*sqrt((1+beta(k,weak))*Hweak);
            a_s = 2*y(k,strong)*sqrt((1+beta(k,strong))*Hstrong);

            d_s = y(k,strong)^2 * Hstrong;

            qx = a_w/(2*sqrt(x)) - a_s/(2*sqrt(1-x)) + d_s;

            fprintf('k=%d q(alpha_weak)=%.4e\n', k, qx);
    end % k loop
 
end
















% function [beta,zeta,y,z] = get_constants_ldt_qt1(...
%     para, H, H_c,decoding_order, alpha,w_k,eta)
%    % Parameters
%     K   = para.K;
%     K_c = para.K_c;
%     M   = para.M;
%     % eta = para.eta;
%     P_max = para.P_max;
%     noise = para.noise;

%     intra_i = zeros(K, K_c); % Initialize A_n vector
%     inteer_i = zeros(K, K_c); % Initialize B_n vector
%     inteer_b = zeros(K, K_c); % Initialize B_n vector
%     inteer_b_all = zeros(K, K_c); % Initialize B_n vector

%     A = zeros(K, K_c); % Initialize A_n vector
%     B = zeros(K, K_c); % Initialize B_n vector
%     A_c = zeros(K,1); % Initialize A_f vector
%     B_c = zeros(K,1); % Initialize B_f vector

%     beta = zeros(K, K_c);
%     y = zeros(K, K_c);
%     zeta = zeros(K, 1);
%     z = zeros(K, 1);

%         %% ===== MAIN LOOP =====
%         for k = 1:K

%             order_k = decoding_order(k,:); % weak → strong (ascending order)
            
%             % Identify strong user (last user in decoding order)
%             strong_user = order_k(end);  % Strong user has highest index in decoding order
            
%             for i = 1:K_c
                
%                 %% ---------- INTER-CLUSTER INTERFERENCE ----------
%                 inter = 0;
%                 inter_b = 0;
%                 h_jk_c=0;

%                 for j = 1:K
%                     if j ~= k
%                         inter   = inter  + abs(H{k,i}*w_k(:, j)).^2;
%                         inter_b = inter_b + abs(H_c{k,i}*w_k(:, j)).^2*eta(j);
%                         h_jk_c = h_jk_c + abs(H_c{k,i}*w_k(:, j)).^2;
%                     end
%                 end

%                 %% ---------- INTRA-CLUSTER INTERFERENCE (NOMA) ----------
%                 intra = 0;
%                 pos = find(order_k == i);
%                 H_i = abs(H{k,i}*w_k(:, k)).^2;

%                 for idx = pos+1:K_c
%                     j_user = order_k(idx);
%                     intra = intra + alpha(k,j_user) * H_i;
%                 end

%                 %% ---------- NOMA SIGNAL ----------
%                 signal =  abs(H{k,i}*w_k(:, k)).^2 * alpha(k,i);
%                 A(k,i)= signal ;

%                 %% ---------- NOMA INTERFERENCE ----------
%                 B(k,i) = itnra + inter + inter_b ...
%                           + abs(H_c{k,i}*w_k(:, k)).^2* eta(k) ...
%                           + noise;

%                 beta(k,i)=A(k,i)/B(k,i);
%                 y(k,i)=sqrt((1+beta(k,i))*A(k,i))/(A(k,i)+B(k,i) );

%                 %% ===== BACKSCATTER CONSTRAINTS (ONLY FOR STRONG USER) =====
%                 if i == strong_user
%                     A_c(k) = abs(H_c{k,i}*w_k(:, k)).^2*eta(k); 

%                         % Backscatter interference constraint
%                     B_c(k)  = inter + inter_b + noise;
%                     zeta(k) = A_c(k)/B_c(k);
%                     z(k)=sqrt((1+zeta(k))*A_c(k))/( A_c(k)+B_c(k));
%                 end
%             end
%         end

% end
