function eta = update_eta(para, H, H_c, decoding_order, ...
                           alpha, w_k, eta, beta, zeta, y, z)
% =========================================================
%  Gauss-Seidel update of eta_k for all k=1..K
%  Uses closed-form solution from LDT+QT derivative:
%
%    eta_k* = clip( z_k^2*(1+zeta_k)*h_{k,i*}^c / C_k^2,
%                   eta_k_min, eta_k_max )
%
%  where C_k = sum_i y_{k,i}^2 * h_{k,i}^c + z_k^2 * h_{k,i*}^c
% =========================================================
 
    K     = para.K;
    K_c   = para.K_c;
    noise = para.noise;
    gamma_min = 2^para.R_min_n - 1;
    gamma_c_min = 2^para.R_min_c - 1;
    
 
    for k = 1:K
 
        order_k     = decoding_order(k,:);
        strong_user = order_k(end);   % i* for cluster k
 
        %% ---- Compute C_k ----
        % C_k = sum_{i=1}^{K_c} y_{k,i}^2 * h_{k,i}^c
        %       + z_k^2 * h_{k,i*}^c
        C_k = 0;
        H_c_kk = 0;
 
        for i = 1:K_c
            H_c_ki = abs(H_c{k,i} * w_k(:,k)).^2;   % h_{k,i}^c
            C_k    = C_k + y(k,i)^2 * H_c_ki;
 
            if i == strong_user
                H_c_kk = H_c_ki;   % h_{k,i*}^c — save for later
            end
        end
        C_k = C_k + z(k)^2 * H_c_kk;
 
        %% ---- Unconstrained optimum ----
        % eta_k_unc = z_k^2 * (1+zeta_k) * h_{k,i*}^c / C_k^2
        eta_unc = z(k)^2 * (1 + zeta(k)) * H_c_kk / C_k^2;

        disp(['Unc Eta :', num2str(eta)]);
 
        %% ---- Lower bound from (C2): BD QoS constraint ----
        % eta_k_min = gamma_kc_min * ( sum_{j~=k} h_{j,i*}^c*eta_j
        %                              + I_{k,i*}^{btot} ) / h_{k,i*}^c
        %
        % Gauss-Seidel: use already-updated eta(j) for j<k,
        %               previous eta(j) for j>k (naturally in eta vector)
        inter_b_others = 0;
        inter_others   = 0;
 
        for j = 1:K
            if j ~= k
                H_c_jk = abs(H_c{k,strong_user} * w_k(:,j)).^2; % h_{j,i*}^c
                inter_b_others = inter_b_others + H_c_jk * eta(j);
                inter_others   = inter_others   + abs(H{k,strong_user} * w_k(:,j)).^2;
            end
        end
 
        I_btot    = inter_others + noise;   % I_{k,i*}^{btot}
        eta_k_min = gamma_c_min * (inter_b_others + I_btot) / H_c_kk;
 
        %% ---- Upper bound from (C1): NOMA QoS constraint ----
        % eta_k_max = min(1, (h_{k,i}*alpha_{k,i} - gamma_min*I_tot)
        %                    / (gamma_min * h_{k,i}^c)
        %                    - sum_{j~=k} h_{j,i}^c/h_{k,i}^c * eta_j )
        %
        % Taken over the most binding subcarrier i
        eta_k_max = 1;   % start from hardware limit
 
        for i = 1:K_c
            H_ki  = abs(H{k,i}   * w_k(:,k)).^2;   % h_{k,i}
            H_c_ki = abs(H_c{k,i} * w_k(:,k)).^2;   % h_{k,i}^c (j=k term)
 
            % Compute intra-cluster interference for user i
            pos   = find(order_k == i);
            intra = 0;
            for idx = pos+1:K_c
                j_user = order_k(idx);
                intra  = intra + alpha(k, j_user) * H_ki;
            end
 
            % I_{k,i}^{tot} = intra + inter (no BD interference)
            inter = 0;
            inter_b_others_noma = 0;
            for j = 1:K
                if j ~= k
                    inter = inter + abs(H{k,i} * w_k(:,j)).^2;
                    inter_b_others_noma = inter_b_others_noma ...
                        + abs(H_c{k,i} * w_k(:,j)).^2 * eta(j);
                end
            end
            I_tot = intra + inter + noise;
 
            % Rearranged (C1) for eta_k:
            % h_{k,i}*alpha_{k,i} >= gamma_min*(sum_k h_{k,i}^c*eta_k + I_tot)
            % => eta_k <= (h_{k,i}*alpha_{k,i}/gamma_min - I_tot
            %              - sum_{j~=k} h_{j,i}^c*eta_j) / h_{k,i}^c
            numerator = H_ki * alpha(k,i) / gamma_min ...
                        - I_tot ...
                        - inter_b_others_noma;
            eta_k_max_i = numerator / H_c_ki;
 
            % Take most binding constraint across all subcarriers
            eta_k_max = min(eta_k_max, eta_k_max_i);
        end
 
        %% ---- Ensure valid interval ----
        eta_k_min = max(0, eta_k_min);          % never below 0
        eta_k_max = min(1, max(0, eta_k_max));  % clip to [0,1]
 
        %% ---- Feasibility check ----
        if eta_k_min > eta_k_max
            warning('Cluster %d: infeasible bounds [%.4f, %.4f]. Using eta_k_min.', ...
                    k, eta_k_min, eta_k_max);
            eta_k_max = eta_k_min;
        end
 
        %% ---- Project unconstrained solution onto [eta_min, eta_max] ----
        eta(k) = min(eta_k_max, max(eta_k_min, eta_unc));
 
    end % k loop
 
end