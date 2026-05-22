function [eta, p_alloc, obj_rc_pac] = solve_eta_AO(para, H, H_c, h_c_eh, ...
                                        decoding_order, p_alloc, w_k, eta)

    max_iter = 20;
    epsilon  = para.epsilon;

    obj_rc_pac  = zeros(max_iter, 1);
    true_wsr    = zeros(max_iter, 1);
    
    % Store best feasible solution found
    best_obj    = -inf;
    best_p      = p_alloc;
    best_eta    = eta;

    for n = 1:max_iter
        fprintf('--- Inner iter %d ---\n', n);

        %% Step 1: Update auxiliary variables at CURRENT (p, eta)
        [beta, zeta, y, z, ~, ~, I_tot] = get_constants_ldt_qt1( ...
            para, H, H_c, decoding_order, p_alloc, w_k, eta);

        %% Step 2: Solve CVX subproblem
        [p_new, eta_new, obj_lbt] = ldt_qt_cvx_pac_rc( ...
            para, H, H_c, h_c_eh, I_tot, beta, y, zeta, z, decoding_order, w_k);

        %% Step 3: Check CVX succeeded before accepting solution
        if obj_lbt <= 0
            warning('CVX failed at iter %d, keeping previous point', n);
            % Do NOT update p_alloc/eta — stay at current point
            obj_rc_pac(n) = obj_rc_pac(max(n-1,1));
            true_wsr(n)   = true_wsr(max(n-1,1));
            continue;   % skip to next iteration
        end

        %% Step 4: Verify SCA monotonicity
        if n > 1 && obj_lbt < obj_rc_pac(n-1) - 1e-4
            warning('SCA objective DECREASED: %.6f -> %.6f', ...
                     obj_rc_pac(n-1), obj_lbt);
            % This indicates auxiliary variable update bug
        end

        %% Accept new point
        p_alloc = p_new;
        eta     = eta_new;

        %% Step 5: Evaluate TRUE rate (for logging only)
        [sum_rate, R, R_c] = compute_wsr( ...
            para, H, H_c, decoding_order, p_alloc, w_k, eta);

        %% Track best feasible point
        if sum_rate > best_obj
            best_obj = sum_rate;
            best_p   = p_alloc;
            best_eta = eta;
        end

        obj_rc_pac(n) = obj_lbt;    % SCA objective for convergence
        true_wsr(n)   = sum_rate;   % True rate for logging

        fprintf('Iter %d | SCA obj = %.6f | True WSR = %.6f\n', ...
                 n, obj_lbt, sum_rate);
        fprintf('SCA vs True gap = %.6f\n', obj_lbt - sum_rate);
        disp('p_alloc:'); disp(p_alloc);
        disp('eta:');     disp(eta');
        disp('R:');       disp(R);
        disp('R_c:');     disp(R_c);

        %% Step 6: Convergence on SCA objective
        if n > 1
            rel_change = abs(obj_rc_pac(n) - obj_rc_pac(n-1)) / ...
                         max(1e-6, abs(obj_rc_pac(n-1)));

            fprintf('Rel change (SCA) = %.2e\n', rel_change);

            if rel_change < epsilon
                fprintf('Converged at iter %d\n', n);
                obj_rc_pac = obj_rc_pac(1:n);
                break;
            end
        end
    end

    %% Return best feasible point found (not just last iterate)
    p_alloc = best_p;
    eta     = best_eta;
    
    fprintf('\nFinal: best true WSR = %.6f\n', best_obj);
end









% function [eta,p_alloc,obj_rc_pac] = solve_eta_AO(para, H, H_c, h_c_eh, decoding_order, p_alloc, w_k, eta)

%     max_iter = 20;
%     epsilon  = para.epsilon;


%     obj_rc_pac = zeros(max_iter, 1);

%     for n = 1:max_iter
%        disp(['=====================================================', num2str(n)]);
%         % Step 1: update LDT/QT auxiliary variables
%         [beta, zeta, y, z, ~, ~, I_tot] = get_constants_ldt_qt1( ...
%             para, H, H_c, decoding_order, p_alloc, w_k, eta);
         

%         % % Step 2: jointly update alpha and eta
%         [p_alloc, eta, obj_lbt] = ldt_qt_cvx_pac_rc( ...
%             para, H, H_c, h_c_eh, I_tot, beta, y, zeta, z, decoding_order, w_k);
%         % gamma_noma = 2^para.R_min_n - 1;

%     %     [p_alloc, obj_lbt] = optimize_alpha_cvx( ...
%     % para, H, H_c, I_tot, beta, y, decoding_order, w_k, eta);

            


%         % Step 3: evaluate true WSR
%         [sum_rate, R, R_c, ~, ~, ~, ~, intra_i, inteer_i, inteer_b, inteer_b_all] = compute_wsr( ...
%             para, H, H_c, decoding_order, p_alloc, w_k, eta);

%         obj_rc_pac(n) = sum_rate;

%         fprintf('\nIteration %d\n', n);
%         fprintf('True WSR = %.6f, LDT/QT objective = %.6f\n', sum_rate, obj_lbt);

%         disp('Power allocation alpha:');
%         disp(p_alloc);

%         disp('Reflection coefficients eta:');
%         disp(eta);

%         disp('NOMA rates R:');
%         disp(R);

%         disp('BD rates R_c:');
%         disp(R_c);

%         if n > 1
%             diff_obj = obj_rc_pac(n) - obj_rc_pac(n-1);
%             rel_change = abs(diff_obj) / max(1, abs(obj_rc_pac(n-1)));

%             fprintf('WSR change = %.6e\n', diff_obj);

%             if rel_change < epsilon
%                 obj_rc_pac = obj_rc_pac(1:n);
%                 break;
%             end
%         end

%         disp(['=====================================================', num2str(n)]);
%     end

% end