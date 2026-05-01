function [p_alloc, obj] = optimize_alpha_cvx( ...
    para, H, H_c, I_tot, beta, t, decoding_order, w_k, eta)
% =========================================================
%  Optimise power allocation p(alpha) only.
%  eta is fixed.
% =========================================================

    K          = para.K;
    K_c        = para.K_c;
    gamma_noma = 2^para.R_min_n - 1;

    %% --- Precompute channel scalars ---
    h_ki    = zeros(K, K_c);
    h_c_tot = zeros(K, K_c, K);

    for k = 1:K
        for i = 1:K_c
            h_ki(k,i) = abs(H{k,i} * w_k(:,k))^2;

            for l = 1:K
                h_c_tot(k,i,l) = abs(H_c{k,i} * w_k(:,l))^2;
            end
        end
    end

    %% ===== CVX PROBLEM: alpha only =====
    cvx_begin quiet
        % cvx_solver mosek
        % cvx_precision high


        variable p(K, K_c) nonnegative

        sum_rate = 0;

        for k = 1:K
            order_k = decoding_order(k,:);

            for i = 1:K_c

                %% BD interference is fixed because eta is fixed
                bd_interf_noma = 0;
                for l = 1:K
                    bd_interf_noma = bd_interf_noma + h_c_tot(k,i,l) * eta(l);
                end

                %% Intra-cluster NOMA interference
                pos = find(order_k == i);
                intra = 0;

                for idx = pos+1:K_c
                    j_user = order_k(idx);
                    intra = intra + p(k,j_user);
                end

                %% A, B, A+B
                A_ki = h_ki(k,i) * p(k,i);
                B_ki = I_tot(k,i) + bd_interf_noma + h_ki(k,i) * intra;
                ApB  = A_ki + B_ki;

                %% NOMA LDT + QT objective only
                sum_rate = sum_rate ...
                    + log2(1 + beta(k,i)) - beta(k,i) ...
                    + 2*t(k,i)*sqrt((1 + beta(k,i))*h_ki(k,i))*sqrt(p(k,i)) ...
                    - t(k,i)^2 * ApB;

                %% NOMA QoS constraint
                h_ki(k,i)*p(k,i) >= gamma_noma * ...
                    (h_ki(k,i)*intra + bd_interf_noma + I_tot(k,i));

            end
        end

        %% PAC sum constraint
        for k = 1:K
            sum(p(k,:)) == 1;
        end

        maximize(sum_rate);

    cvx_end

    p_alloc = double(p);
    obj     = cvx_optval;

    fprintf('Alpha-only CVX status: %s | obj = %.6f\n', cvx_status, obj);

    if ~strcmp(cvx_status,'Solved') && ~strcmp(cvx_status,'Inaccurate/Solved')
        warning('Alpha-only CVX did not solve cleanly: %s', cvx_status);
    end
end