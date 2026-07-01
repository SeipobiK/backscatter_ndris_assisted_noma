function [Rates,obj_history, w_k_history, theta_history, ...
          rate_f_history, rate_n_history, rate_c_history, ...
          noma_signal_history, BST_signal_history, ...
          noma_interference_history, BST_interference_history, ...
          intra_cluster_history, inter_cluster_history, ...
          inter_cluster_BST_history, inter_cluster_BST_all_history, ...
          decoding_order_history, alpha_history,gains_it_history,eta_history] = bd_ris_opt(para, channel_data, J_r, J_t)

    
    outer_iter = para.outer_iter;
    N = para.N;
    M = para.M;
    K = para.K;
    K_c = para.K_c;
    delta = zeros(N, N);
    grad_prev = zeros(N, N);
    
    % Initialize storage
    obj_history = zeros(outer_iter+1, 1);
    w_k_history = zeros(M, K, outer_iter+1);
    theta_history = zeros(N, N, outer_iter+1);
    rate_f_history = zeros(K, outer_iter+1);
    rate_n_history = zeros(K, outer_iter+1);
    rate_c_history = zeros(K, outer_iter+1);
    eta_history = zeros(K, outer_iter+1);
    
    % Interference and signal 
    noma_signal_history = zeros(K, K_c, outer_iter+1);
    Rates = zeros(K, K_c, outer_iter+1);
    BST_signal_history = zeros(K, outer_iter+1);
    noma_interference_history = zeros(K, K_c, outer_iter+1);
    BST_interference_history = zeros(K, outer_iter+1);
    intra_cluster_history = zeros(K, K_c, outer_iter+1);
    inter_cluster_history = zeros(K, K_c, outer_iter+1);
    inter_cluster_BST_history = zeros(K, K_c, outer_iter+1);
    inter_cluster_BST_all_history = zeros(K, K_c, outer_iter+1);
    decoding_order_history = zeros(K, K_c, outer_iter+1);
    alpha_history = zeros(K, K_c, outer_iter+1);
    gains_it_history = zeros(K, K_c, outer_iter+1);
    h_c_eh = cell(K, K_c);

    
    % Initialize RIS phase shifts
    % Theta = initialize_ris_phases(para);

    Theta = initialize_bd_ris(para, 'random');
    % Theta = diag(Theta);
    
    % Build effective channels
    [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);
    
    % Initialize beamforming vectors
    W_init = initialize_beamforming(para, H);
    % w_k=initialize_beamforming(para, H);
    
    % Determine decoding order
    [decoding_order,gains_it] = determine_decoding_order(para, H, W_init);
        % Update decoding order
    % decoding_order = determine_decoding_order(para, H, w_k);
    alpha = assign_power_by_decoding_order(para, decoding_order);
    
    % Initialize Taylor parameters
    [A, B, Ac, Bc] = initialize_taylor_parameters(para);

    eta=para.eta;
    
    %% Initial feasibility optimization
    % Active BF feasibility
    [W_opt, A_abf, B_abf, Ac_abf, Bc_abf, obj_prev, status] = ...
        optimize_feasibility_abf(para, H, H_c, channel_data, decoding_order, A, B, Ac, Bc, alpha,eta);

    w_k = extract_beamforming_vectors(W_opt);



    
    % Update effective channels
    [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);

    [decoding_order,gains_it] = determine_decoding_order(para, H, w_k);
    alpha = assign_power_by_decoding_order(para, decoding_order);
    gains_it_history(:, :, 1) = gains_it;
    alpha_history(:, :, 1) = alpha;

    % Compute initial metrics (iteration 0)
    [sum_rate, R, R_c, noma_signal, BST_signal, noma_interference, BST_interference, ...
        intra_i, inteer_i, inteer_b, inteer_b_all] = compute_wsr(...
        para, H, H_c, decoding_order,alpha, w_k,eta);
    
    % Store initial results
    obj_history(1) = sum_rate;
    rate_f_history(:, 1) = R(:, 2);
    rate_n_history(:, 1) = R(:, 1);
    rate_c_history(:, 1) = R_c;
    eta_history(:, 1) = para.eta(1);
    w_k_history(:, :, 1) = w_k;
    theta_history(:, :, 1) = Theta;
    
    % Store interference and signal metrics
    Rates(:, :, 1) = R;
    noma_signal_history(:, :, 1) = noma_signal;
    BST_signal_history(:, 1) = BST_signal;
    noma_interference_history(:, :, 1) = noma_interference;
    BST_interference_history(:, 1) = BST_interference;
    intra_cluster_history(:, :, 1) = intra_i;
    inter_cluster_history(:, :, 1) = inteer_i;
    inter_cluster_BST_history(:, :, 1) = inteer_b;
    inter_cluster_BST_all_history(:, :, 1) = inteer_b_all;
    decoding_order_history(:, :, 1) = decoding_order;
    
    % disp(['Initial WSR = ', num2str(sum_rate)]);

    
    %% Main iteration loop
    for iter = 1:outer_iter
          
        [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t); 

        %%  Active beamforming optimization-Algorithm  1
        [W_opt, A_abf, B_abf, Ac_abf, Bc_abf, obj_curr, status, ~] = ...
            maximize_sum_rate_iterative_abf(para, H, H_c, channel_data, decoding_order, ...
            A_abf, B_abf, Ac_abf, Bc_abf,alpha,eta);       
        
        w_k = extract_beamforming_vectors(W_opt); 
   
            % Update effective channels
            [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);

            % [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);

            for k = 1:K
                for i = 1:K_c                
                    % Backscatter link: H_c = g_b' * J_t' * Theta' * J_r' * H_all * f
                    h_c_eh{k,i} = channel_data.g_b{k}' * J_r * Theta * J_t * ...
                                channel_data.H_all ;
                end
            end
            % disp(['After ABF and PBF optimization, before eta update: WSR = ', num2str(sum_rate_)]);   

        % %% ===== INITIALIZE DUAL VARIABLES =====
        % % NOMA users (K clusters × K_c users per cluster)
        % duals_vars_init = zeros(para.K, para.K_c);
        % % Alternative (if you want faster initial constraint enforcement):
        % % duals_vars_init = 0.01 * ones(para.K, para.K_c);

        % % Backscatter (one per cluster)
        % duals_vars_init_c = zeros(para.K, 1);

        % % Learning rate (fixed for now)
        % learning_rate = 0.01;

        % %% ===== CALL YOUR FUNCTION =====
        % [beta, zeta, y, z, sum_rate_const_alpha, sum_rate_const_eta, ...
        % duals_vars_init, duals_vars_init_c] = ...
        %     compute_fp_gual_vars(para, H, H_c, decoding_order, alpha, w_k, eta, ...
        %     duals_vars_init, learning_rate, duals_vars_init_c)

        % [grad] = compute_eucl_grad(para,w_k,Theta,channel_data,decoding_order,...
        %    beta,zeta,y,z,duals_vars_init,duals_vars_init_c,alpha,eta);


           [sum_rate, R, R_c, noma_signal, BST_signal, noma_interference, BST_interference, ...
            intra_i, inteer_i, inteer_b, inteer_b_all] = compute_wsr(...
            para, H, H_c, decoding_order, alpha, w_k,eta);


            disp(sum_rate);
            % disp(Theta'*Theta);


[Theta_opt, nu_opt, duals_vars, duals_vars_c] = ...
    passive_beamforming_riemannian(para, w_k, Theta, channel_data, ...
                                   decoding_order, alpha, eta);
       [H, H_c] = build_effective_channels(para, channel_data, Theta_opt, J_r, J_t);  

               %% Compute metrics
        [sum_rate, R, R_c, noma_signal, BST_signal, noma_interference, BST_interference, ...
            intra_i, inteer_i, inteer_b, inteer_b_all] = compute_wsr(...
            para, H, H_c, decoding_order, alpha, w_k,eta);

        disp(['After ABF, PBF and passive optimization, WSR = ', num2str(sum_rate)]);
        % disp(Theta_opt'*Theta_opt );
fff
        
        %% Update decoding order
        [decoding_order,gains_it]= determine_decoding_order(para, H, w_k);
        gains_it_history(:, :, iter+1) = gains_it;

        alpha_history(:, :, iter+1) = alpha;

        

        jifdhjfhuiosdhfgj

        % Store resultsss
        obj_history(iter+1) = sum_rate;
        rate_f_history(:, iter+1) = R(:, 2);
        rate_n_history(:, iter+1) = R(:, 1);
        rate_c_history(:, iter+1) = R_c;
        w_k_history(:, :, iter+1) = w_k;
        theta_history(:, :, iter+1) = Theta;
        
        % Store interference and signal metrics
        noma_signal_history(:, :, iter+1) = noma_signal;
        BST_signal_history(:, iter+1) = BST_signal;
        noma_interference_history(:, :, iter+1) = noma_interference;
        BST_interference_history(:, iter+1) = BST_interference;
        intra_cluster_history(:, :, iter+1) = intra_i;
        inter_cluster_history(:, :, iter+1) = inteer_i;
        inter_cluster_BST_history(:, :, iter+1) = inteer_b;
        inter_cluster_BST_all_history(:, :, iter+1) = inteer_b_all;
        decoding_order_history(:, :, iter+1) = decoding_order;
        Rates(:, :, iter+1) = R;
        eta_history(:, iter+1) = eta;

    end
    
    % disp(['Size of obj_history: ', num2str(size(obj_history))]);
end
