function [eta,p_alloc,obj_rc_pac] = solve_eta_AO(para, H, H_c, h_c_eh, decoding_order, p_alloc, w_k, eta)

    max_iter = 20;
    epsilon  = para.epsilon;


    obj_rc_pac = zeros(max_iter, 1);

    for n = 1:max_iter
       disp(['=====================================================', num2str(n)]);
        % Step 1: update LDT/QT auxiliary variables
        [beta, zeta, y, z, ~, ~, I_tot] = get_constants_ldt_qt1( ...
            para, H, H_c, decoding_order, p_alloc, w_k, eta);
         

        % % Step 2: jointly update alpha and eta
        [p_alloc, eta, obj_lbt] = ldt_qt_cvx_pac_rc( ...
            para, H, H_c, h_c_eh, I_tot, beta, y, zeta, z, decoding_order, w_k);
        % gamma_noma = 2^para.R_min_n - 1;

    %     [p_alloc, obj_lbt] = optimize_alpha_cvx( ...
    % para, H, H_c, I_tot, beta, y, decoding_order, w_k, eta);

            


        % Step 3: evaluate true WSR
        [sum_rate, R, R_c, ~, ~, ~, ~, intra_i, inteer_i, inteer_b, inteer_b_all] = compute_wsr( ...
            para, H, H_c, decoding_order, p_alloc, w_k, eta);

        obj_rc_pac(n) = sum_rate;

        fprintf('\nIteration %d\n', n);
        fprintf('True WSR = %.6f, LDT/QT objective = %.6f\n', sum_rate, obj_lbt);

        disp('Power allocation alpha:');
        disp(p_alloc);

        disp('Reflection coefficients eta:');
        disp(eta);

        disp('NOMA rates R:');
        disp(R);

        disp('BD rates R_c:');
        disp(R_c);

        if n > 1
            diff_obj = obj_rc_pac(n) - obj_rc_pac(n-1);
            rel_change = abs(diff_obj) / max(1, abs(obj_rc_pac(n-1)));

            fprintf('WSR change = %.6e\n', diff_obj);

            if rel_change < epsilon
                obj_rc_pac = obj_rc_pac(1:n);
                break;
            end
        end

        disp(['=====================================================', num2str(n)]);
    end

end