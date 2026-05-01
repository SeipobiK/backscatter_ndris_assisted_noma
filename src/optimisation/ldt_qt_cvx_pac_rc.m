function [p_alloc, eta, obj] = ldt_qt_cvx_pac_rc(...
    para, H, H_c, h_c_eh, I_tot, beta, t, zeta, s, decoding_order, w_k)
% =========================================================
%  Joint optimisation of p (alpha) and eta via CVX/MOSEK.
%
%  IMPORTANT: I_tot(k,i) passed in must equal inter+noise only
%             (as returned by get_constants_ldt_qt1).
%             Do NOT add inter_bf again inside this function.
% =========================================================

    K          = para.K;
    K_c        = para.K_c;
    gamma_noma = 2^para.R_min_n - 1;
    gamma_back = 2^para.R_min_c - 1;
    eh_min     = para.bst_threshold;
    rho_eh     = para.rho;

    %% --- Precompute ALL channel scalars outside CVX ---
    % This avoids any DCP parsing issues inside cvx_begin
    h_ki    = zeros(K, K_c);   % h_{k,i}
    h_c_kik = zeros(K, K_c);   % h_{k,i}^c  own BD k at user i

    % h_{k,i}^{c,tot}_l = |H_c{k,i}*w_k(:,l)|^2  for each (k,i,l)
    h_c_tot = zeros(K, K_c, K);  % h_c_tot(k,i,l) = BD l interf at user (k,i)

    % h_{k,i*}^{c,int}_l = |H_c{k,i*}*w_k(:,l)|^2 for l~=k (BD interf at BD k)
    h_c_int = zeros(K, K);       % h_c_int(k,l) = BD l interf at BD k strong user

    for k = 1:K
        order_k     = decoding_order(k,:);
        strong_user = order_k(end);
        for i = 1:K_c
            h_ki(k,i)    = abs(H{k,i}   * w_k(:,k))^2;
            h_c_kik(k,i) = abs(H_c{k,i} * w_k(:,k))^2;
            for l = 1:K
                h_c_tot(k,i,l) = abs(H_c{k,i} * w_k(:,l))^2;
            end
            if i == strong_user
                for l = 1:K
                    if l ~= k
                        h_c_int(k,l) = abs(H_c{k,i} * w_k(:,l))^2;
                    end
                end
            end
        end
    end

    %% ===== CVX PROBLEM =====
    cvx_begin quiet
        cvx_solver mosek


        variable p(K, K_c) nonnegative
        variable eta_var(K, 1) nonnegative

        sum_rate = 0;

        for k = 1:K
            order_k     = decoding_order(k,:);
            strong_user = order_k(end);

            for i = 1:K_c

                %% --- BD interference at NOMA user (k,i) from ALL BDs ---
                % = sum_{l=1}^{K} h_c_tot(k,i,l) * eta(l)
                % I_tot(k,i) already contains inter_bf + noise
                % so total denominator = I_tot(k,i) + bd_interf + h_ki*intra
                bd_interf_noma = 0;
                for l = 1:K
                    bd_interf_noma = bd_interf_noma ...
                        + h_c_tot(k,i,l) * eta_var(l);
                end

                %% --- Intra-cluster NOMA interference ---
                pos   = find(order_k == i);
                intra = 0;
                for idx = pos+1:K_c
                    j_user = order_k(idx);
                    intra  = intra + p(k, j_user);
                end

                %% --- Signal and interference terms ---
                % A_{k,i} = h_{k,i} * p(k,i)
                % B_{k,i} = I_tot(k,i) + bd_interf + h_{k,i}*intra
                % A+B     = I_tot(k,i) + bd_interf + h_{k,i}*(intra+p(k,i))
                A_ki = h_ki(k,i) * p(k,i);
                B_ki = I_tot(k,i) + bd_interf_noma + h_ki(k,i) * intra;
                ApB  = A_ki + B_ki;

                %% --- NOMA objective term ---
                sum_rate = sum_rate ...
                    + log2(1 + beta(k,i)) - beta(k,i) ...
                    + 2 * t(k,i) * sqrt((1+beta(k,i)) * h_ki(k,i)) * sqrt(p(k,i)) ...
                    - t(k,i)^2  * ApB;

                %% --- NOMA QoS (C1) ---
                % h_{k,i}*p(k,i) >= gamma*(h_{k,i}*intra + bd_interf + I_tot)
                h_ki(k,i)*p(k,i) >= gamma_noma * ...
                    (h_ki(k,i)*intra + bd_interf_noma + I_tot(k,i));

                %% --- BD terms: strong user only ---
                if i == strong_user

                    %% BD interference from OTHER BDs at strong user
                    % = sum_{l~=k} h_c_int(k,l) * eta(l)
                    % I_tot(k,i) already has inter_bf + noise for this user
                    bd_interf_BD = 0;
                    for l = 1:K
                        if l ~= k
                            bd_interf_BD = bd_interf_BD ...
                                + h_c_int(k,l) * eta_var(l);
                        end
                    end

                    %% BD signal and denominator
                    % A_k^c = h_{k,i*}^c * eta(k)
                    % B_k^c = bd_interf_BD + I_tot(k,i*)
                    % A+B   = h_{k,i*}^c*eta(k) + bd_interf_BD + I_tot
                    A_kc  = h_c_kik(k,i) * eta_var(k);
                    B_kc  = bd_interf_BD + I_tot(k,i);
                    ApB_c = A_kc + B_kc;

                    %% BD objective term
                    sum_rate = sum_rate ...
                        + log2(1+zeta(k)) - zeta(k) ...
                        + 2*s(k)*sqrt((1+zeta(k))*h_c_kik(k,i))*sqrt(eta_var(k))...
                        - s(k)^2 * ApB_c;

                    %% BD QoS (C2)
                    h_c_kik(k,i)*eta_var(k) >= gamma_back * ...
                        (bd_interf_BD + I_tot(k,i));

                    %% Energy harvesting constraint
                    if eh_min > 0
                        h_eh = abs(h_c_eh{k,i} * w_k(:,k))^2;
                        (1 - eta_var(k)) * rho_eh * h_eh >= eh_min;
                    end

                end
            end % i
        end % k

        %% Power budget (C3)
        for k = 1:K
            sum(p(k,:)) == 1;
                    %% Box constraint (C4)
            0 <= eta_var(k) <= 1;
        end



        maximize(sum_rate);

    cvx_end

    p_alloc = double(p);
    eta     = double(eta_var);
    obj     = cvx_optval;

    fprintf('CVX status: %s | obj = %.6f\n', cvx_status, obj);

    if ~strcmp(cvx_status,'Solved') && ~strcmp(cvx_status,'Inaccurate/Solved')
        warning('CVX did not solve cleanly: %s', cvx_status);
    end
end