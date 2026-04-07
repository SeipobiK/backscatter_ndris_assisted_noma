function [V_opt, A_opt, B_opt, Ac_opt, Bc_opt, obj_curr, status] = optimize_feasibility_pbf(para,w_k,channel_data,decoding_order,...
     A, B, Ac, Bc,J_r,J_t,alpha)
    % Perform feasibility optimization with iterative updates
    
    max_iterations = para.max_iter;
        Ac_opt = Ac;
    Bc_opt = Bc;
    A_opt = A;
    B_opt = B;

    obj_history = zeros(max_iterations, 1);

    
    for iter = 1:max_iterations
        % Call feasibility optimization solver
        % Note: You need to ensure your 'feasible' function accepts these parameters
   [V_opt, A_opt, B_opt, Ac_opt, Bc_opt, obj_curr, status] =sca_feasible_pbf(para,w_k,channel_data,decoding_order,...
    A_opt, B_opt, Ac_opt, Bc_opt,alpha,J_r,J_t);

    % fprintf('Iter %3d | Objective: %.4e', iter, obj_curr);
        
        % Check solver status
        if ~strcmp(status, 'Solved')
            warning('Feasibility optimization failed at iteration %d, status: %s', iter, status);
            break;
        end
        
        % Update history and display progress
        obj_history(iter) = obj_curr;
        % display_optimization_progress(iter, obj_curr);
        % disp(['Iteration pbf  ',num2str(iter) ]);
        % disp(['obj_curr pbf  ',num2str(obj_curr) ]);
        
        % Check convergence
        if iter > 1 && abs(obj_history(iter)) < 1e-9
            obj_history = obj_history(1:iter);
            % fprintf('  Feasibility optimization converged after %d iterations\n', iter);
            % disp(['obj_curr  ',num2str(obj_history(iter)) ]);
            break;
        end
    end
end
