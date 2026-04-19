function [alpha_opt, eta_opt, history] = max_sumrate_pac_rc(para, H_eff, H_back, ...
    I_total, gamma_min_noma, gamma_min_back, decoding_order, max_iter, tol)
    % fp_optimization - Main FP algorithm for joint optimization
    % Input:
    %   para - simulation parameters
    %   H_eff - effective channel gains h_{k,i}
    %   H_back - backscatter channel gains h_{k,i}^c
    %   I_total - total interference + noise
    %   gamma_min_noma - minimum SINR for NOMA users (K x K_c)
    %   gamma_min_back - minimum SINR for backscatter (K x 1)
    %   decoding_order - SIC decoding order (K x K_c)
    %   max_iter - maximum number of iterations
    %   tol - convergence tolerance
    % Output:
    %   alpha_opt - optimal power allocation (K x K_c)
    %   eta_opt - optimal reflection coefficients (K x 1)
    %   history - structure with convergence history
    
    K = para.K;
    K_c = para.K_c;
    
    % Initialize feasible alpha and eta
    alpha = zeros(K, K_c);
    for k = 1:K
        % Equal power allocation initially
        alpha(k, :) = 1/K_c * ones(1, K_c);
    end
    eta = 0.5 * ones(K, 1);  % Start with 0.5 reflection
    
    % Initialize auxiliary variables
    beta = zeros(K, K_c);
    t = zeros(K, K_c);
    zeta = zeros(K, 1);
    s = zeros(K, 1);
    
    % Store history
    history.obj = zeros(max_iter, 1);
    history.alpha = cell(max_iter, 1);
    history.eta = cell(max_iter, 1);
    history.beta = cell(max_iter, 1);
    history.t = cell(max_iter, 1);
    history.zeta = cell(max_iter, 1);
    history.s = cell(max_iter, 1);
    
    % Main FP loop
    for iter = 1:max_iter
        % Store previous values for convergence check
        alpha_prev = alpha;
        eta_prev = eta;
        
        % Step 1: Update auxiliary variables
        [beta, t, zeta, s] = update_auxiliary_variables(para, H_eff, H_back, ...
            I_total, alpha, eta, decoding_order);
        
        % Step 2: Solve main convex subproblem
        [alpha, eta, status] = solve_cvx_subproblem(para, H_eff, H_back, I_total, ...
            beta, t, zeta, s, gamma_min_noma, gamma_min_back, decoding_order);
        
        % Compute objective value for convergence monitoring
        obj_value = compute_objective(para, H_eff, H_back, I_total, ...
            alpha, eta, beta, t, zeta, s, decoding_order);
        
        % Store history
        history.obj(iter) = obj_value;
        history.alpha{iter} = alpha;
        history.eta{iter} = eta;
        history.beta{iter} = beta;
        history.t{iter} = t;
        history.zeta{iter} = zeta;
        history.s{iter} = s;
        
        % Check convergence
        if iter > 1
            alpha_change = norm(alpha(:) - alpha_prev(:), 'inf');
            eta_change = norm(eta - eta_prev, 'inf');
            obj_change = abs(history.obj(iter) - history.obj(iter-1));
            
            fprintf('Iter %d: obj=%.6f, alpha_change=%.2e, eta_change=%.2e\n', ...
                iter, obj_value, alpha_change, eta_change);
            
            if max(alpha_change, eta_change) < tol && obj_change < tol
                fprintf('Converged at iteration %d\n', iter);
                history.converged = true;
                history.iterations = iter;
                break;
            end
        end
        
        % Check CVX status
        if ~strcmp(status, 'Solved')
            fprintf('Warning: CVX status = %s at iteration %d\n', status, iter);
        end
    end
    
    % Set outputs
    alpha_opt = alpha;
    eta_opt = eta;
    
    % Trim history
    history.obj = history.obj(1:iter);
    history.alpha = history.alpha(1:iter);
    history.eta = history.eta(1:iter);
    history.beta = history.beta(1:iter);
    history.t = history.t(1:iter);
    history.zeta = history.zeta(1:iter);
    history.s = history.s(1:iter);
end

% Helper function to compute objective
function obj = compute_objective(para, H_eff, H_back, I_total, alpha, eta, beta, t, zeta, s, decoding_order)
    K = para.K;
    K_c = para.K_c;
    
    % NOMA rates
    noma_obj = 0;
    for k = 1:K
        for i = 1:K_c
            intra_interf = 0;
            for j = 1:K_c
                if decoding_order(k, j) >= decoding_order(k, i)
                    intra_interf = intra_interf + alpha(k, decoding_order(k, j));
                end
            end
            
            denom = I_total(k,i) + H_back(k,i) * sum(eta) + H_eff(k,i) * intra_interf;
            
            noma_obj = noma_obj + log2(1+beta(k,i)) - beta(k,i) ...
                + 2*t(k,i) * sqrt((1+beta(k,i)) * H_eff(k,i) * alpha(k,i)) ...
                - t(k,i)^2 * denom;
        end
    end
    
    % Backscatter rates
    back_obj = 0;
    for k = 1:K
        i_star = decoding_order(k, end);
        denom = I_total(k, i_star) + H_back(k, i_star) * sum(eta);
        
        back_obj = back_obj + log2(1+zeta(k)) - zeta(k) ...
            + 2*s(k) * sqrt((1+zeta(k)) * H_back(k, i_star) * eta(k)) ...
            - s(k)^2 * denom;
    end
    
    obj = noma_obj + back_obj;
end