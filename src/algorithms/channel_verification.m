function [Rates,obj_history, w_k_history, theta_history, ...
          rate_f_history, rate_n_history, rate_c_history, ...
          noma_signal_history, BST_signal_history, ...
          noma_interference_history, BST_interference_history, ...
          intra_cluster_history, inter_cluster_history, ...
          inter_cluster_BST_history, inter_cluster_BST_all_history, ...
          decoding_order_history, alpha_f, alpha_n] = channel_verification(para, channel_data, J_r, J_t)
    
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
        A_abf, B_abf, Ac_abf, Bc_abf, J_r, J_t, alpha);
    
    % Extract RIS phase shifts
    Theta = extract_theta(V_opt, para, H, H_c, decoding_order, alpha, w_k);
    
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
    
    % disp(['Initial WSR = ', num2str(sum_rate)]);

    h=zeros(K, K_c);
    h_c=zeros(K, K_c);
    h_c_eh=zeros(K, K_c);
    
    %% Main iteration loop
    for iter = 1:outer_iter
        eta=para.eta;

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
    %    verify_intra_interference(decoding_order, noma_signal, alpha, intra_i);
                

        [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);


        % Compute WSR
        [sum_rate, R, R_c, noma_signal, BST_signal, noma_interference, BST_interference, ...
            intra_i, I_tot, inteer_b, inteer_b_all] = compute_wsr(...
            para, H, H_c, decoding_order, alpha, w_k);



                %%  Active beamforming optimization-Algorithm  1
        [W_opt, A_abf, B_abf, Ac_abf, Bc_abf, obj_curr, status, ~] = ...
            maximize_sum_rate_iterative_abf(para, H, H_c, channel_data, decoding_order, ...
            A_abf, B_abf, Ac_abf, Bc_abf,alpha);       
        
        w_k = extract_beamforming_vectors(W_opt); 

                %% Passive beamforming optimization-Algorithm 2
        [V_opt, A_pbf, B_pbf, Ac_pbf, Bc_pbf, ~, SR, converged] = ...
            maximize_sum_rate_iterative_pbf(para, w_k, channel_data, decoding_order, ...
            A_pbf, B_pbf, Ac_pbf, Bc_pbf,alpha, J_r, J_t);
        
            % % Extract RIS phase shifts
            Theta = extract_theta(V_opt, para, H, H_c, decoding_order,alpha, w_k);
   
            % Update effective channels
            [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);



        
        % for i=1:10
        %      disp(['sum rate before :', num2str(sum_rate)]);
        %     % Update auxiliary varables
        %     [beta,I_tot,t,h,h_c,h_eh] = update_beta_zeta_st(...
        %         para, H, H_c, decoding_order, alpha, w_k, eta);

        %     % ldt/qt pac optimisation
        %     [alpha, cvx_status, obj] = ldt_find_eta(...
        %         para, h,R_c,I_tot, beta, t, decoding_order);
        %     % sum_rate = compute_sum_rate(para, H, H_c, decoding_order, alpha, eta, w_k);

        %      % Compute WSR
        %      [sum_rate, R, R_c, noma_signal, BST_signal, noma_interference, BST_interference, ...
        %     intra_i, I_tot, inteer_b, inteer_b_all] = compute_wsr(...
        %     para, H, H_c, decoding_order, alpha, w_k);

        %     for j=1:2     
        %         disp(['beta :',num2str(beta(j,1)), ' alpha :',num2str(alpha(2,1)), ' SINR from compute_wsr :', num2str(inteer_b(j,1)),' obj form sum_wsr :', num2str(sum_rate)]);
        %         disp([' beta :', num2str(beta(j,2))]);
        %     end
        % end


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
        
    end
    
    for c=1:K
        alpha_f(c) = para.alpha_k_f;
        alpha_n(c) = para.alpha_k_n; 
    end
    % disp(['Size of obj_history: ', num2str(size(obj_history))]);
end






% function [Rates,obj_history, w_k_history, theta_history, ...
%           rate_f_history, rate_n_history, rate_c_history, ...
%           noma_signal_history, BST_signal_history, ...
%           noma_interference_history, BST_interference_history, ...
%           intra_cluster_history, inter_cluster_history, ...
%           inter_cluster_BST_history, inter_cluster_BST_all_history, ...
%           decoding_order_history, alpha_f, alpha_n] = channel_verification(para, channel_data, J_r, J_t)
    
%     outer_iter = para.outer_iter;
%     N = para.N;
%     M = para.M;
%     K = para.K;
%     K_c = para.K_c;
    
%     % Initialize storage
%     obj_history = zeros(outer_iter+1, 1);
%     w_k_history = zeros(M, K, outer_iter+1);
%     theta_history = zeros(N, N, outer_iter+1);
%     rate_f_history = zeros(K, outer_iter+1);
%     rate_n_history = zeros(K, outer_iter+1);
%     rate_c_history = zeros(K, outer_iter+1);
    
%     % Interference and signal 
%     noma_signal_history = zeros(K, K_c, outer_iter+1);
%     Rates = zeros(K, K_c, outer_iter+1);
%     BST_signal_history = zeros(K, outer_iter+1);
%     noma_interference_history = zeros(K, K_c, outer_iter+1);
%     BST_interference_history = zeros(K, outer_iter+1);
%     intra_cluster_history = zeros(K, K_c, outer_iter+1);
%     inter_cluster_history = zeros(K, K_c, outer_iter+1);
%     inter_cluster_BST_history = zeros(K, K_c, outer_iter+1);
%     inter_cluster_BST_all_history = zeros(K, K_c, outer_iter+1);
%     decoding_order_history = zeros(K, K_c, outer_iter+1);
    
%     % Initialize RIS phase shifts
%     Theta = initialize_ris_phases(para);
%     % Theta = diag(Theta);
    
%     % Build effective channels
%     [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);
    
%     % Initialize beamforming vectors
%     W_init = initialize_beamforming(para, H);
%     w_k=initialize_beamforming(para, H);
    
%     % Determine decoding order
%     decoding_order = determine_decoding_order(para, H, W_init);
%         % Update decoding order
%     % decoding_order = determine_decoding_order(para, H, w_k);
%     alpha = assign_power_by_decoding_order(para, decoding_order);
    
%     % Initialize Taylor parameters
%     [A, B, Ac, Bc] = initialize_taylor_parameters(para);
    
%     %% Initial feasibility optimization
%     % Active BF feasibility
%     [W_opt, A_abf, B_abf, Ac_abf, Bc_abf, obj_prev, status] = ...
%         optimize_feasibility_abf(para, H, H_c, channel_data, decoding_order, A, B, Ac, Bc, alpha);
    
%     w_k = extract_beamforming_vectors(W_opt);
    
%     % Passive BF feasibility
%     [V_opt, A_pbf, B_pbf, Ac_pbf, Bc_pbf, obj_curr, status] = ...
%         optimize_feasibility_pbf(para, w_k, channel_data, decoding_order, ...
%         A, B, Ac, Bc, J_r, J_t, alpha);
    
%     % Extract RIS phase shifts
%     Theta = extract_theta(V_opt, para, H, H_c, decoding_order, alpha, w_k);
%     Theta = eye(N); 
    
%     % Update effective channels
%     [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);
    
%     % Update decoding order
%     decoding_order = determine_decoding_order(para, H, w_k);
%     alpha = assign_power_by_decoding_order(para, decoding_order);
    
%     % Compute initial metrics (iteration 0)
%     [sum_rate, R, R_c, noma_signal, BST_signal, noma_interference, BST_interference, ...
%         intra_i, inteer_i, inteer_b, inteer_b_all] = compute_wsr(...
%         para, H, H_c, decoding_order,alpha, w_k);
    
%     % Store initial results
%     obj_history(1) = sum_rate;
%     rate_f_history(:, 1) = R(:, 2);
%     rate_n_history(:, 1) = R(:, 1);
%     rate_c_history(:, 1) = R_c;
%     w_k_history(:, :, 1) = w_k;
%     theta_history(:, :, 1) = Theta;
    
%     % Store interference and signal metrics
%     Rates(:, :, 1) = R;
%     noma_signal_history(:, :, 1) = noma_signal;
%     BST_signal_history(:, 1) = BST_signal;
%     noma_interference_history(:, :, 1) = noma_interference;
%     BST_interference_history(:, 1) = BST_interference;
%     intra_cluster_history(:, :, 1) = intra_i;
%     inter_cluster_history(:, :, 1) = inteer_i;
%     inter_cluster_BST_history(:, :, 1) = inteer_b;
%     inter_cluster_BST_all_history(:, :, 1) = inteer_b_all;
%     decoding_order_history(:, :, 1) = decoding_order;
    
%     disp(['Initial WSR = ', num2str(sum_rate)]);
    
%     %% Main iteration loop
%     for iter = 1:outer_iter
%         % Active beamforming optimization
%         disp(['active beamforming Vectors... ', num2str(iter)]);
%         disp(w_k);
%         disp('active beamforming Vectors...');


%         disp('Passive beamforming Vectors...');
%         disp(diag(Theta));
%         disp('Passive beamforming Vectors...');


%         disp(['Ennd BF... ', num2str(iter)]);


%         disp(['Iteration======================= ', num2str(iter), ': Status = ', status]);
%         for k=1:K
%             disp(['Cluster ', num2str(k)]);
%             for i=1:K_c  
%                 sinr=2^R(k, i)-1;

%                  disp(['User  ', num2str(i),' Alpha = ', num2str(alpha(k, i)), ' : noma_signal = ', num2str(noma_signal(k, i)) ...
%                  , ' : Intra interference ', num2str(intra_i(k, i)), ' User rate = : ', num2str(R(k, i)),  ': SINR =  ', num2str(sinr)]);      
%             end   
%        end
%        disp('Decoding order: ');
%        disp(decoding_order);
%        disp(['Iteration======================= ', num2str(iter), ': Status = ', status]);   
%        verify_intra_interference(decoding_order, noma_signal, alpha, intra_i);
%        alpha = assign_power_by_decoding_order(para, decoding_order);
         
%         [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);

%         %  Active beamforming optimization-Algorithm  1
%         [W_opt, A_abf, B_abf, Ac_abf, Bc_abf, obj_curr, status, ~] = ...
%             maximize_sum_rate_iterative_abf(para, H, H_c, channel_data, decoding_order, ...
%             A_abf, B_abf, Ac_abf, Bc_abf,alpha);
%         disp(['Iteration ', num2str(iter), ': Status = ', status]);
       
        
%         w_k = extract_beamforming_vectors(W_opt);
%         % disp(w_k);
        
%         % % Passive beamforming optimization-Algorithm 2
%         [V_opt, A_pbf, B_pbf, Ac_pbf, Bc_pbf, ~, ~, converged] = ...
%             maximize_sum_rate_iterative_pbf(para, w_k, channel_data, decoding_order, ...
%             A_pbf, B_pbf, Ac_pbf, Bc_pbf,alpha, J_r, J_t);
        
%         % Extract RIS phase shifts
%         Theta = extract_theta(V_opt, para, H, H_c, decoding_order,alpha, w_k);
   
%         % Update effective channels
%         [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);
        
%         % Update decoding order
%         % decoding_order = determine_decoding_order(para, H, w_k);

        
%         % Compute metrics
%         [sum_rate, R, R_c, noma_signal, BST_signal, noma_interference, BST_interference, ...
%             intra_i, inteer_i, inteer_b, inteer_b_all] = compute_wsr(...
%             para, H, H_c, decoding_order, alpha, w_k);

%         % disp(['Iteration ', num2str(iter), ': noma_signal = ', num2str(sum_rate)]);
%         % disp(noma_signal);

%         decoding_order = determine_decoding_order(para, H, w_k);


    
%         % Store resultsss
%         obj_history(iter+1) = sum_rate;
%         rate_f_history(:, iter+1) = R(:, 2);
%         rate_n_history(:, iter+1) = R(:, 1);
%         rate_c_history(:, iter+1) = R_c;
%         w_k_history(:, :, iter+1) = w_k;
%         theta_history(:, :, iter+1) = Theta;
        
%         % Store interference and signal metrics
%         noma_signal_history(:, :, iter+1) = noma_signal;
%         BST_signal_history(:, iter+1) = BST_signal;
%         noma_interference_history(:, :, iter+1) = noma_interference;
%         BST_interference_history(:, iter+1) = BST_interference;
%         intra_cluster_history(:, :, iter+1) = intra_i;
%         inter_cluster_history(:, :, iter+1) = inteer_i;
%         inter_cluster_BST_history(:, :, iter+1) = inteer_b;
%         inter_cluster_BST_all_history(:, :, iter+1) = inteer_b_all;
%         decoding_order_history(:, :, iter+1) = decoding_order;
%         Rates(:, :, iter+1) = R;
        
%         disp(['Outer Iteration ', num2str(iter), ': WSR = ', num2str(sum_rate)]);
        
%         % % Check convergence (optional)
%         % if iter > 1 && abs(obj_history(iter+1) - obj_history(iter)) < 1e-4
%         %     disp(['Converged at iteration ', num2str(iter)]);
%         %     % Truncate history to actual iterations
%         %     obj_history = obj_history(1:iter+1);
%         %     w_k_history = w_k_history(:, :, 1:iter+1);
%         %     theta_history = theta_history(:, :, 1:iter+1);
%         %     rate_f_history = rate_f_history(:, 1:iter+1);
%         %     rate_n_history = rate_n_history(:, 1:iter+1);
%         %     rate_c_history = rate_c_history(:, 1:iter+1);
%         %     noma_signal_history = noma_signal_history(:, :, 1:iter+1);
%         %     BST_signal_history = BST_signal_history(:, 1:iter+1);
%         %     noma_interference_history = noma_interference_history(:, :, 1:iter+1);
%         %     BST_interference_history = BST_interference_history(:, 1:iter+1);
%         %     intra_cluster_history = intra_cluster_history(:, :, 1:iter+1);
%         %     inter_cluster_history = inter_cluster_history(:, :, 1:iter+1);
%         %     inter_cluster_BST_history = inter_cluster_BST_history(:, :, 1:iter+1);
%         %     inter_cluster_BST_all_history = inter_cluster_BST_all_history(:, :, 1:iter+1);
%         %     decoding_order_history = decoding_order_history(:, :, 1:iter+1);
%         %     Rates(:, :, iter+1) = R;
%         %     break;
%         % end
%     end
    
%     for c=1:K
%         alpha_f(c) = para.alpha_k_f;
%         alpha_n(c) = para.alpha_k_n; 
%     end
%     % disp(['Size of obj_history: ', num2str(size(obj_history))]);
% end