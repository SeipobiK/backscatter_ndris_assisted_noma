function [V_opt, A_opt, B_opt, A_c_opt, B_c_opt, obj_history, obj_new, converged] = ...
    maximize_sum_rate_iterative_pbf(para, w_k, channel_data, decoding_order, ...
    A_opt, B_opt, A_c_opt, B_c_opt, alpha, J_r, J_t)

    %% Parameters
    max_iter   = 25;
    N          = para.N;
    MC_MAX     = para.MC_MAX;

    %% Initialization
    obj_history    = zeros(max_iter, 1);
    obj_history_mc = zeros(max_iter, MC_MAX);
    converged      = false;

    step_size      = zeros(max_iter, 1);
    relax_parameter= zeros(max_iter, 1);

    max_eig_val    = zeros(max_iter, 1);
    max_eig_vec    = zeros(N, max_iter);
    U_opt          = zeros(N, N, max_iter);

    %% ---- Initial SCA solve ----
    [V_opt, A_opt, B_opt, A_c_opt, B_c_opt, obj_prev, status] = ...
        sca_rate_max_pbf_init(para, w_k, channel_data, decoding_order, ...
        A_opt, B_opt, A_c_opt, B_c_opt, alpha, J_r, J_t);

    if ~strcmp(status, 'Solved')
        warning('Initial SCA problem not solved');
        obj_history = NaN;
        return;
    end

    %% Eigen decomposition (iteration 1)
    [v_max, lambda_max] = max_eigVect(V_opt);

    U_opt(:,:,1)     = V_opt;
    max_eig_val(1)   = lambda_max;
    max_eig_vec(:,1) = v_max;
    obj_history(1)   = obj_prev;

    %% Initial relaxation
    ratio_0           = lambda_max / trace(V_opt);
    step_size(1)      = 0.45 * (1 - ratio_0);
    relax_parameter(1)= min(1, ratio_0 + step_size(1));

    %% Display initial status
    % fprintf('Iter 1 | Obj: %.6f | λ_max: %.4e | Ratio: %.4f | Rank≈%d\n', ...
    %         obj_history(1), lambda_max, ratio_0, rank(V_opt));

    %% ---- Main Iteration Loop ----
    for m = 2:max_iter

        %% Solve SCA subproblem
        [V_new, A_new, B_new, A_c_new, B_c_new, obj_new, cvx_status] = ...
            sca_rate_max_pbf(para, w_k, channel_data, decoding_order, ...
            A_opt, B_opt, A_c_opt, B_c_opt, alpha, J_r, J_t, ...
            max_eig_vec(:,m-1), relax_parameter(m-1));

        if strcmp(cvx_status, 'Solved')

            %% Accept update
            V_opt = V_new;
            A_opt = A_new;
            B_opt = B_new;
            A_c_opt = A_c_new;
            B_c_opt = B_c_new;

            obj_history(m) = obj_new;

            %% Eigen update
            [v_max, eig_max] = max_eigVect(V_opt);

            U_opt(:,:,m)     = V_opt;
            max_eig_val(m)   = eig_max;
            max_eig_vec(:,m) = v_max;

            %% Full eigenvalue spectrum
            eig_vals = sort(real(eig(V_opt)), 'descend');
            dominant_eig = eig_vals(1);
            eig_ratio    = dominant_eig / sum(eig_vals);
            eig_gap      = dominant_eig - eig_vals(min(2, length(eig_vals)));

            %% Step size update
            step_size(m) = 0.45* (1 - eig_ratio);

            % %% Display progress every iteration
            fprintf('Iter %2d | Obj: %.6f | ΔObj: %.3e | λ_max: %.4e | Ratio: %.4f | Gap: %.4e | Rank≈%d\n', ...
                    m, obj_history(m), obj_history(m)-obj_history(m-1), ...
                    dominant_eig, eig_ratio, eig_gap, rank(V_opt));
            % disp(eig_vals');

        else
            %% Failure handling
            U_opt(:,:,m)   = U_opt(:,:,m-1);
            max_eig_val(m) = max_eig_val(m-1);
            obj_history(m) = obj_history(m-1);

            step_size(m) = step_size(m-1) / 2;

            warning('SCA failed at iteration %d, reducing step size', m);

            if step_size(m) < 1e-6
                obj_history = NaN;
                return;
            end

            continue;
        end

        %% Relaxation parameter update
        ratio = max_eig_val(m) / trace(U_opt(:,:,m));
        relax_parameter(m) = min(1, ratio + 0.45 * (1 - ratio));

        %% Convergence check
        if abs(obj_history(m) - obj_history(m-1)) < 1e-3 && ...
           abs(1 - relax_parameter(m)) <= 1e-5

            converged   = true;
            obj_history = obj_history(1:m);
            fprintf('Converged at iteration %d\n', m);
            break;
        end
    end
end