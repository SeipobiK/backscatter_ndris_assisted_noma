function [W_opt, A_opt, B_opt, Ac_opt, Bc_opt, obj_curr, status,obj_history] = maximize_sum_rate_iterative_abf(para,H, H_c,channel_data,decoding_order,...
     A, B, Ac, Bc,alpha)
    % Perform feasibility optimization with iterative updates
    
    max_iterations = para.max_iter;

    A_opt = A;
    B_opt = B;

    obj_history = zeros(max_iterations, 1);

    
    Ac_opt = Ac;
    Bc_opt = Bc;
    for iter = 1:max_iterations
        % Call feasibility optimization solver
        % Note: You need to ensure your 'feasible' function accepts these parameters
   [W_opt, A_opt, B_opt, Ac_opt, Bc_opt, obj_curr, status] = sca_rate_max_abf(para,channel_data, H, H_c, A_opt, B_opt, Ac_opt, Bc_opt, decoding_order, alpha);
        
        % Check solver status
        if ~strcmp(status, 'Solved')
            warning('Active BF  at iteration %d, status: %s', iter, status);
            break;
        end
        
        % Update history and display progress
        obj_history(iter) = obj_curr;
        % display_optimization_progress(iter, obj_curr);
        % if iter > 1
        %     disp(['Iteration abf ',num2str(iter), ' Difference: ', num2str(abs(obj_history(iter)-obj_history(iter-1))) ] );
        %     disp(['obj_curr abf ',num2str(obj_curr) ]);     
        % end

        % Check convergence
        if iter > 1 && abs(obj_history(iter)-obj_history(iter-1)) < 1e-3
            obj_history = obj_history(1:iter);
            % fprintf('  Active BF solved after %d iterations\n', iter);
            % disp(['obj_curr  ',num2str(obj_history(iter))]);
            break;
        end
    end
end
