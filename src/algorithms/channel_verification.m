function [Rates,obj_history, w_k_history, theta_history, ...
          rate_f_history, rate_n_history, rate_c_history, ...
          noma_signal_history, BST_signal_history, ...
          noma_interference_history, BST_interference_history, ...
          intra_cluster_history, inter_cluster_history, ...
          inter_cluster_BST_history, inter_cluster_BST_all_history, ...
          decoding_order_history, alpha_f, alpha_n] = channel_verification(para, channel_data, J_r, J_t)
    % run_optimization - Main optimization function for both DRIS and NDRIS
    % Input:
    %   para - simulation parameters
    %   channel_data - preprocessed channel data
    %   J_r - receiving RIS matrix (identity for DRIS, optimized for NDRIS)
    %   J_t - transmitting RIS matrix (identity for DRIS, optimized for NDRIS)
    % Output:
    %   obj_history - WSR history (outer_iter+1 x 1)
    %   w_k_history - beamforming vectors history (M x K x (outer_iter+1))
    %   theta_history - RIS phase shifts history (N x N x (outer_iter+1))
    %   rate_f_history - far user rates history (K x (outer_iter+1))
    %   rate_n_history - near user rates history (K x (outer_iter+1))
    %   rate_c_history - common rate history (K x (outer_iter+1))
    %   noma_signal_history - NOMA signal power history (K x K_c x (outer_iter+1))
    %   BST_signal_history - BST signal power history (K x (outer_iter+1))
    %   noma_interference_history - NOMA interference history (K x K_c x (outer_iter+1))
    %   BST_interference_history - BST interference history (K x (outer_iter+1))
    %   intra_cluster_history - intra-cluster interference (K x K_c x (outer_iter+1))
    %   inter_cluster_history - inter-cluster interference (K x K_c x (outer_iter+1))
    %   inter_cluster_BST_history - inter-cluster BST interference (K x K_c x (outer_iter+1))
    %   inter_cluster_BST_all_history - total inter-cluster BST interference (K x K_c x (outer_iter+1))
    %   decoding_order_history - decoding order history (K x K_c x (outer_iter+1))
    %   alpha_f - far user power allocation
    %   alpha_n - near user power allocation
    
    outer_iter = para.outer_iter;
    N = para.N;
    M = para.M;
    K = para.K;
    K_c = para.K_c;
    
    % Initialize storage
    obj_history = zeros(outer_iter+1, 1);
    w_k_history = zeros(M, K, outer_iter+1);
    theta_history = zeros(N, N, outer_iter+1);
    rate_f_history = zeros(K, outer_iter+1);
    rate_n_history = zeros(K, outer_iter+1);
    rate_c_history = zeros(K, outer_iter+1);
    
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
    
    % Initialize RIS phase shifts
    Theta = initialize_ris_phases(para);
    % Theta = diag(Theta);
    
    % Build effective channels
    [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);
    
    % Initialize beamforming vectors
    W_init = initialize_beamforming(para, H);
    w_k=initialize_beamforming(para, H);
    
    % Determine decoding order
    decoding_order = determine_decoding_order(para, H, W_init);
        % Update decoding order
    % decoding_order = determine_decoding_order(para, H, w_k);
    alpha = assign_power_by_decoding_order(para, decoding_order);
    
    % Initialize Taylor parameters
    [A, B, Ac, Bc] = initialize_taylor_parameters(para);
    
    %% Initial feasibility optimization
    % Active BF feasibility
    [W_opt, A_abf, B_abf, Ac_abf, Bc_abf, obj_prev, status] = ...
        optimize_feasibility_abf(para, H, H_c, channel_data, decoding_order, A, B, Ac, Bc, alpha);
    
    w_k = extract_beamforming_vectors(W_opt);
    
    % Passive BF feasibility
    [V_opt, A_pbf, B_pbf, Ac_pbf, Bc_pbf, obj_curr, status] = ...
        optimize_feasibility_pbf(para, w_k, channel_data, decoding_order, ...
        A, B, Ac, Bc, J_r, J_t, alpha);
    
    % Extract RIS phase shifts
    Theta = extract_theta(V_opt, para, H, H_c, decoding_order, alpha, w_k);
    Theta = eye(N); 
    
    % Update effective channels
    [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);
    
    % Update decoding order
    decoding_order = determine_decoding_order(para, H, w_k);
    alpha = assign_power_by_decoding_order(para, decoding_order);
    
    % Compute initial metrics (iteration 0)
    [sum_rate, R, R_c, noma_signal, BST_signal, noma_interference, BST_interference, ...
        intra_i, inteer_i, inteer_b, inteer_b_all] = compute_wsr(...
        para, H, H_c, decoding_order,alpha, w_k);
    
    % Store initial results
    obj_history(1) = sum_rate;
    rate_f_history(:, 1) = R(:, 2);
    rate_n_history(:, 1) = R(:, 1);
    rate_c_history(:, 1) = R_c;
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
    
    disp(['Initial WSR = ', num2str(sum_rate)]);
    
    %% Main iteration loop
    for iter = 1:outer_iter
        % Active beamforming optimization
        disp(['active beamforming Vectors... ', num2str(iter)]);
        disp(w_k);
        disp('active beamforming Vectors...');


        disp('Passive beamforming Vectors...');
        disp(diag(Theta));
        disp('Passive beamforming Vectors...');


        disp(['Ennd BF... ', num2str(iter)]);


        disp(['Iteration======================= ', num2str(iter), ': Status = ', status]);
        for k=1:K
            disp(['Cluster ', num2str(k)]);
            for i=1:K_c  
                sinr=2^R(k, i)-1;

                 disp(['User  ', num2str(i),' Alpha = ', num2str(alpha(k, i)), ' : noma_signal = ', num2str(noma_signal(k, i)) ...
                 , ' : Intra interference ', num2str(intra_i(k, i)), ' User rate = : ', num2str(R(k, i)),  ': SINR =  ', num2str(sinr)]);      
            end   
       end
       disp('Decoding order: ');
       disp(decoding_order);
       disp(['Iteration======================= ', num2str(iter), ': Status = ', status]);   
       verify_intra_interference(decoding_order, noma_signal, alpha, intra_i);
       alpha = assign_power_by_decoding_order(para, decoding_order);
         
  err

        [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);

        %  Active beamforming optimization-Algorithm  1
        [W_opt, A_abf, B_abf, Ac_abf, Bc_abf, obj_curr, status, ~] = ...
            maximize_sum_rate_iterative_abf(para, H, H_c, channel_data, decoding_order, ...
            A_abf, B_abf, Ac_abf, Bc_abf,alpha);
        disp(['Iteration ', num2str(iter), ': Status = ', status]);
       
        
        w_k = extract_beamforming_vectors(W_opt);
        % disp(w_k);
        
        % % Passive beamforming optimization-Algorithm 2
        [V_opt, A_pbf, B_pbf, Ac_pbf, Bc_pbf, ~, ~, converged] = ...
            maximize_sum_rate_iterative_pbf(para, w_k, channel_data, decoding_order, ...
            A_pbf, B_pbf, Ac_pbf, Bc_pbf,alpha, J_r, J_t);
        
        % Extract RIS phase shifts
        Theta = extract_theta(V_opt, para, H, H_c, decoding_order,alpha, w_k);
   
        % Update effective channels
        [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);
        
        % Update decoding order
        % decoding_order = determine_decoding_order(para, H, w_k);

        
        % Compute metrics
        [sum_rate, R, R_c, noma_signal, BST_signal, noma_interference, BST_interference, ...
            intra_i, inteer_i, inteer_b, inteer_b_all] = compute_wsr(...
            para, H, H_c, decoding_order, alpha, w_k);

        % disp(['Iteration ', num2str(iter), ': noma_signal = ', num2str(sum_rate)]);
        % disp(noma_signal);

        decoding_order = determine_decoding_order(para, H, w_k);


    
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
        
        disp(['Outer Iteration ', num2str(iter), ': WSR = ', num2str(sum_rate)]);
        
        % % Check convergence (optional)
        % if iter > 1 && abs(obj_history(iter+1) - obj_history(iter)) < 1e-4
        %     disp(['Converged at iteration ', num2str(iter)]);
        %     % Truncate history to actual iterations
        %     obj_history = obj_history(1:iter+1);
        %     w_k_history = w_k_history(:, :, 1:iter+1);
        %     theta_history = theta_history(:, :, 1:iter+1);
        %     rate_f_history = rate_f_history(:, 1:iter+1);
        %     rate_n_history = rate_n_history(:, 1:iter+1);
        %     rate_c_history = rate_c_history(:, 1:iter+1);
        %     noma_signal_history = noma_signal_history(:, :, 1:iter+1);
        %     BST_signal_history = BST_signal_history(:, 1:iter+1);
        %     noma_interference_history = noma_interference_history(:, :, 1:iter+1);
        %     BST_interference_history = BST_interference_history(:, 1:iter+1);
        %     intra_cluster_history = intra_cluster_history(:, :, 1:iter+1);
        %     inter_cluster_history = inter_cluster_history(:, :, 1:iter+1);
        %     inter_cluster_BST_history = inter_cluster_BST_history(:, :, 1:iter+1);
        %     inter_cluster_BST_all_history = inter_cluster_BST_all_history(:, :, 1:iter+1);
        %     decoding_order_history = decoding_order_history(:, :, 1:iter+1);
        %     Rates(:, :, iter+1) = R;
        %     break;
        % end
    end
    
    for c=1:K
        alpha_f(c) = para.alpha_k_f;
        alpha_n(c) = para.alpha_k_n; 
    end
    % disp(['Size of obj_history: ', num2str(size(obj_history))]);
end