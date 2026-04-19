function [p_alloc, cvx_status, obj] = ldt_find_eta(...
    para, H,R_c,I_tot, beta, t, decoding_order)

    % Parameters
    K = para.K;
    K_c = para.K_c;
    R_min_noma = 2^para.R_min_n - 1;      % SINR target for NOMA users

    cvx_begin quiet
        % cvx_solver mosek
        
        % Variables (using p for power allocation to avoid conflicts)
        variable p(K, K_c) nonnegative
        
        % Objective
        expression sum_rate
        sum_rate = 0;
        for k = 1:K
            order_k = decoding_order(k, :);  % weak → strong
            for i = 1:K_c

                % Intra-cluster interference (from stronger users)
                pos = find(order_k == i);
                intra = 0;
                for idx = pos+1:K_c
                    j_user = order_k(idx);
                    intra = intra + p(k, j_user);
                end
                
                % Objective: NOMA user rate (using LDT + QT)
                sum_rate = sum_rate + log2(1 + beta(k,i)) - beta(k,i) ...
                    + 2 * t(k,i) * sqrt((1 + beta(k,i)) * H(k,i) * p(k,i)) ...
                    - t(k,i)^2 * (I_tot(k,i) + H(k,i) * (intra + p(k,i))) ;

                % disp(2 * t(k,i) * sqrt((1 + beta(k,i)) * H(k,i) * p(k,i)) );
                % disp(t(k,i)^2 * (I_tot(k,i) + H(k,i) * (intra + p(k,i))) );


                % QoS constraint for NOMA user
                H(k,i) * p(k,i) >= R_min_noma * (I_tot(k,i) + H(k,i) * intra);

            end
            sum_rate=sum_rate+R_c(k);
        end
        
        % Power allocation constraints (sum of p per cluster = 1)
        for k = 1:K
            sum(p(k, :)) == 1;
        end
        
        % Maximize objective
        maximize(sum_rate)
        
    cvx_end
    
    % Output
    obj = cvx_optval;
    p_alloc = p;  % Output renamed to p_alloc
    
    % Check feasibility
    if ~strcmp(cvx_status, 'Solved')
        fprintf('Warning: CVX status = %s\n', cvx_status);
    end
end