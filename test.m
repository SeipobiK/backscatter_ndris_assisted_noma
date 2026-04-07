load("results/2026/mar/30/full_workspace_dris_M8_N32_20260330_232420.mat");

% disp(decoding_order_dris(1,:,:,1));
% disp(noma_signal_dris(1,:,:,1));

% disp(intra_cluster_dris(1,:,:,1));

disp(rates_dris(1,:,:,1));


% disp(inter_cluster_BST_all_dris(1,:,:,1));






% intra_cluster_dris

% clear all; clc;
% addpath(genpath(pwd));
% rng(2022); % For reproducibility');

% % Initialize parameters
% para = para_init();
% [BS_array, RIS_array] = generate_arrays(para);

% % Constants
% max_feasible = 25;
% max_iter = 15;
% outer_iter = para.outer_iter;
% MC_MAX = para.MC_MAX;
% K = para.K;
% K_c=para.K_c;
% N = para.N;
% M = para.M;
% num_mc_iterations=1;
% K_c=para.K_c;


% % Preallocate results for DRIS and NDRIS only
% obj_history_dris = zeros(outer_iter, MC_MAX);
% obj_history_ndris = zeros(outer_iter, MC_MAX);
% alpha_f_mc_ndris = zeros(K,MC_MAX);
% alpha_n_mc_ndris = zeros(K,MC_MAX);
% w_k_ndris = zeros(M, K,MC_MAX);
% w_k_dris = zeros(M, K,MC_MAX);

% alpha_f_mc_dris = zeros(K,MC_MAX);
% alpha_n_mc_dris = zeros(K,MC_MAX);
 
% rate_f_mc_ndris_outer = zeros(K,outer_iter,MC_MAX);
% rate_n_mc_ndris_outer = zeros(K,outer_iter,MC_MAX);
% rate_c_mc_dris_outer = zeros(K,outer_iter,MC_MAX);


% channel_far_dris= zeros(K,outer_iter,MC_MAX);
% channel_near_dris= zeros(K,outer_iter,MC_MAX);

% inter_cluster_interference_near_dris= zeros(K,outer_iter,MC_MAX);
% inter_cluster_interference_near_bst_dris= zeros(K,outer_iter,MC_MAX);
% inter_cluster_interference_far_dris= zeros(K,outer_iter,MC_MAX);
% inter_cluster_interference_near_b_dris= zeros(K,outer_iter,MC_MAX);

% inter_cluster_interference_dris= zeros(K,K_c,outer_iter,MC_MAX);
% intra_cluster_interference_dris= zeros(K,K_c,outer_iter,MC_MAX);
% inter_cluster_BST_dris= zeros(K,K_c,outer_iter,MC_MAX);
% inter_cluster_BST_all_dris= zeros(K,K_c,outer_iter,MC_MAX);

% signal_dris= zeros(K,K_c,outer_iter,MC_MAX);
% signal_c_dris= zeros(K,outer_iter,MC_MAX);

% total_interference_dris= zeros(K,K_c,outer_iter,MC_MAX);
% total_BST_interference_dris= zeros(K,outer_iter,MC_MAX);

% channel_far_ndris= zeros(K,outer_iter,MC_MAX);
% channel_near_ndris= zeros(K,outer_iter,MC_MAX);

% inter_cluster_interference_near_ndris= zeros(K,outer_iter,MC_MAX);
% inter_cluster_interference_near_bst_ndris= zeros(K,outer_iter,MC_MAX);
% inter_cluster_interference_far_ndris= zeros(K,outer_iter,MC_MAX);
% inter_cluster_interference_near_b_ndris= zeros(K,outer_iter,MC_MAX);

% rate_f_mc_dris_outer = zeros(K,outer_iter,MC_MAX);
% rate_n_mc_dris_outer = zeros(K,outer_iter,MC_MAX);
% rate_f_mc_ndris = zeros(K,MC_MAX);
% rate_n_mc_ndris = zeros(K,MC_MAX);
% rate_c_mc_ndris = zeros(K,MC_MAX);

% rate_f_mc_dris = zeros(K,MC_MAX);
% rate_n_mc_dris = zeros(K,MC_MAX);
% rate_c_mc_dris = zeros(K,MC_MAX);


% rng_seeds = randi(1e6, MC_MAX, 1);

% % Create results directory
% results_dir = 'results_passive';
% if ~exist(results_dir, 'dir')
%     mkdir(results_dir);
% end

% % Start parallel pool
% if isempty(gcp('nocreate'))
%     num_workers = 14;
%     pool = parpool('local', num_workers);
%     fprintf('Using %d workers for parallel processing\n', pool.NumWorkers);
% else
%     pool = gcp;
%     fprintf('Existing pool with %d workers found\n', pool.NumWorkers);
% end

% % Main parallel loop
% tic;
% [BS_array_par, RIS_array_par] = generate_arrays(para);
% %% MAIN MONTE CARLO LOOP
% %% MAIN MONTE CARLO LOOP
% fprintf('Starting Monte Carlo iterations...\n');

% for mc = 1:num_mc_iterations
%     fprintf('Monte Carlo Iteration %d/%d\n', mc, num_mc_iterations);
    
%     % Set random seed for reproducibility
%     rng(rng_seeds(mc), 'twister');
    
%     % Generate channels
%     [H_local, g_local, f_local] = generate_channel(para, BS_array_par, RIS_array_par);
    
%     % Preprocess channels
%     channel_data = preprocess_channels(para, H_local, g_local, f_local);
%     % disp(channel_data.g_b{1});
    
%     % Initialize RIS phase shifts
%     Theta = initialize_ris_phases(para);

%     % % NDRIS: Optimized J matrices
%     % g_LOS_reshaped = reshape(g_local, para.N, []);
%     % J_r = design_J_r(g_LOS_reshaped);
%     % J_t = design_J_t(H_local);

%     J_r = eye(N); % Override with identity for DRIS
%     J_t = eye(N); % Override with identity for DRIS
%     % Determine decoding order
    
%     % Build effective channels
%     [H, H_c] = build_effective_channels(para, channel_data, Theta,J_r,J_t);
    
%     % Initialize beamforming vectors
%     W_init = initialize_beamforming(para, H);
    
%     % J_r = design_J_r(g_LOS_reshaped);
%     % J_t = design_J_t(channel_data.H_all);


%     decoding_order = determine_decoding_order(para, H, W_init);


%     %% Active BF feasible points
%     [A, B, Ac, Bc] = initialize_taylor_parameters(para);
%     [W_opt, A_abf, B_abf, Ac_abf, Bc_abf, obj_prev, status] = optimize_feasibility_abf(para,H, H_c,channel_data,decoding_order,...
%         A, B, Ac, Bc);
    
%       w_k = extract_beamforming_vectors(W_opt);

%     %% Passive BF feasible points
%     [V_opt, A_pbf, B_pbf, Ac_pbf, Bc_pbf, obj_curr, status] = optimize_feasibility_pbf(para,w_k,channel_data,decoding_order,...
%      A_abf, B_abf, Ac_abf, Bc_abf,J_r,J_t);
            
%      Theta = extract_theta( ...
%     V_opt, para, H, H_c,decoding_order,channel_data.alpha,w_k);



%     for iter=1:para.outer_iter

%         %% ACTIVE BEAMFORMING OPTIMIZATION
%         [H, H_c] = build_effective_channels(para, channel_data, Theta,J_r,J_t);

%         [W_opt, A_abf, B_abf, Ac_abf, Bc_abf, obj_curr, status,obj_history] = maximize_sum_rate_iterative_abf(para,H, H_c,channel_data,decoding_order,...
%              A_abf, B_abf, Ac_abf, Bc_abf);
%             %  Extract beamforming vectors
%             w_k = extract_beamforming_vectors(W_opt);
%         %% PASSIVE BEAMFORMING OPTIMIZATION
%         [V_opt, A_pbf, B_pbf, Ac_pbf, Bc_pbf, obj_history_pbf, obj_history_mc, converged] = ...
%         maximize_sum_rate_iterative_pbf(para, w_k, channel_data, decoding_order, ...
%         A_pbf, B_pbf, Ac_pbf, Bc_pbf, channel_data.alpha, J_r, J_t);

%             %  Extract RIS phase shifts
%             [H, H_c] = build_effective_channels(para, channel_data, Theta,J_r,J_t);
%             [Theta, ~] = extract_theta( ...
%                 V_opt, para, H, H_c,decoding_order,channel_data.alpha,w_k);
%             %  Update effective channels
%             [H, H_c] = build_effective_channels(para, channel_data, Theta,J_r,J_t);

%             % update decoding order
%             % decoding_order = determine_decoding_order(para, H, w_k);
%             % Compute WSR and display progress

%             [sum_rate,R,R_c,noma_signal,BST_signal,noma_interference,BST_interference,intra_i,inteer_i,inteer_b,inteer_b_all]  = compute_wsr(...
%                 para, H, H_c,decoding_order, channel_data.alpha,w_k);

%             % [sum_rate,R,R_c,A,A_c,B,B_c,intra_i,inteer_i,inteer_b,inteer_b_all]

%             disp(['Outer Iteration ', num2str(iter), ': WSR = ', num2str(sum_rate)]);
%             obj_history_dris(iter, mc) = sum_rate;

%             % rate_n_mc_dris_outer(:, :, mc) = near_history_dris;
%             % (K,outer_iter,MC_MAX);

%             rate_f_mc_dris_outer(:, iter, mc) = R(:,2);
%             rate_n_mc_dris_outer(:, iter, mc) = R(:,1);
%             rate_c_mc_dris_outer(:, iter, mc) = R_c;

%             channel_far_dris(:, :, mc) = noma_signal(:,2);
%             channel_near_dris(:, :, mc) =noma_signal(:,1);


%             inter_cluster_interference_dris(:,:,iter, mc) = inteer_i;
%             intra_cluster_interference_dris(:,:,iter, mc) = intra_i;
%             inter_cluster_BST_dris(:,:,iter, mc) = inteer_b;
%             inter_cluster_BST_all_dris(:,:,iter, mc) = inteer_b_all;

%             signal_dris(:,:,iter, mc) = noma_signal;
%             signal_c_dris(:,iter, mc) = BST_signal;

%             total_interference_dris(:,:,iter, mc) = noma_interference;
%             total_BST_interference_dris(:,iter, mc) = BST_interference;


%             % inter_cluster_interference_near_dris(:, :, mc) = inter_cluster_interference_n;
%             % inter_cluster_interference_near_bst_dris(:, :, mc) = inter_cluster_interference_n_bst;
%             % inter_cluster_interference_far_dris(:, :, mc) = inter_cluster_interference_f;
%             % inter_cluster_interference_near_b_dris(:, :, mc) = inter_cluster_interference_n_b;

%             disp(['Outer Iteration ', num2str(iter), ': WSR = ', num2str(sum_rate)]);


       
%         % Print summary for this MC run
%         % fprintf('MC %d: Final WSR - DRIS=%.4f, NDRIS=%.4f\n', ...
%         %         mc, obj_history_dris(iter, mc), obj_history_ndris_local(end));

%     end


% end

% disp(obj_history_dris(:, 1))





%     % [sum_rate,R,R_c,~,~,~,~,intra_i,inteer_i,inteer_b,inteer_b_all] = compute_wsr(...
%     % para, H, H_c,decoding_order, channel_data.alpha,w_k);
    

%     % [W_opt, A, B, Ac, Bc, obj_prev, status] = optimize_feasibility_abf(para,H, H_c,channel_data,decoding_order,...
%     %     A, B, Ac, Bc);
    
%     % [W_opt, A_opt, B_opt, Ac_opt, Bc_opt, obj_curr, status,obj_history] = maximize_sum_rate_iterative_abf(para,H, H_c,channel_data,decoding_order,...
%     %  A, B, Ac, Bc);
%     %  disp(obj_history);
%     %  disp(['Final WSR: ', num2str(obj_history(end))]);

%     % w_k = extract_beamforming_vectors(W_opt);



%     % [V_opt, A_opt, B_opt, Ac_opt, Bc_opt, obj_curr, status] = optimize_feasibility_pbf(para,w_k,channel_data,decoding_order,...
%     %  A, B, Ac, Bc,J_r,J_t);

%     % [V_opt, A_opt, B_opt, A_c_opt, B_c_opt, obj_history_pbf, obj_history_mc, converged] = ...
%     % maximize_sum_rate_iterative_pbf(para, w_k, channel_data, decoding_order, ...
%     % A_opt, B_opt, Ac_opt, Bc_opt, channel_data.alpha, J_r, J_t);

%     % [Theta, ~] = extract_theta( ...
%     % V_opt, para, H, H_c,decoding_order,channel_data.alpha,w_k);

%     % [H, H_c] = build_effective_channels(para, channel_data, Theta,J_r,J_t);




%     % [sum_rate,R,R_c,~,~,~,~,intra_i,inteer_i,inteer_b,inteer_b_all] = compute_wsr(...
%     % para, H, H_c,decoding_order, channel_data.alpha,w_k);

%     % % disp(['Final WSR:.Pass ', num2str(obj_history_pbf)]);
%     % disp(['Final WSR:obj_history ', num2str(obj_history(end))]);
%     % disp(['Final WSR:sum_rate ', num2str(sum_rate)]);










    

%     % W = extract_beamforming_vectors(W_opt);
%     % % decoding_order = determine_decoding_order(para, H, W);

%     % % [A, B, Ac, Bc] = initialize_taylor_parameters(para);
%     % disp("computing WSR()444444444444444444444444444444444444");
%     % [W_opt,A_, B_, Ac_, Bc_, obj_prev, status] = sca_rate_max_abf(...
%     % para, H, H_c,  A, B, Ac, Bc,decoding_order, channel_data.alpha);
%     % disp(['Initial WSR AFTER OPTIMIZATION: ', num2str(obj_prev)]);
%     % w_k = extract_beamforming_vectors(W_opt);
%     % % display rank of W_opt
%     % % disp(['Rank of W_opt: ', num2str(rank(W_opt(:,:,1)))]);

%     % [H, H_c] = build_effective_channels(para, channel_data, Theta,J_r,J_t);
%     % disp(['Initial WSR: ', num2str(obj_prev)]);

%     % disp(channel_data.alpha)

%     % disp("computing WSR");


%     % % [A, B, Ac, Bc] = initialize_taylor_parameters(para);
%     % [W_opt,A_, B_, Ac_, Bc_, obj_prev, status] = sca_rate_max_abf(...
%     % para, H, H_c,  A_, B_, Ac_, Bc_,decoding_order, channel_data.alpha);
%     % disp(['WSR after ABF optimization:obj_prev ', num2str(obj_prev)]);
%     % w_k = extract_beamforming_vectors(W_opt);

%     %     [sum_rate,R,R_c,A,A_c,B,B_c,intra_i,inteer_i,inteer_b,inteer_b_all] = compute_wsr(...
%     % para, H, H_c,decoding_order, channel_data.alpha,w_k);
%     % disp(['Initial WSR AFTER OPTIMIZATION:sum_rate ', num2str(sum_rate)]);

%     % % display eigenvalues of W_opt
%     % disp('Eigenvalues of W_opt:');
%     % disp(eig(W_opt(:,:,1)));
%     % disp(eig(W_opt(:,:,2)));

%     % disp(intra_i);

%     % for k=1:K
%     %         disp(['Near user: ',num2str(intra_i(k,1))]);
%     %         disp(['Far user: ',num2str(intra_i(k,2))]);
%     % end
%     % disp(['Rank of W_opt: ', num2str(rank(W_opt))]);
%         % display intra for user that does SIC
%     % for k=1:K
%     %     order_k = decoding_order(k,:); % weak → strong (ascending order)
%     %     i = order_k(end); % Index of the user that performs SIC (strong user)

%     %     disp(['User ', num2str(i), ' performs SIC on user ', num2str(i)]);
%     %     disp(['Intra-cluster interference for user ', num2str(i), ': ', num2str(intra_i(k, i))]);
%     %     disp(['Inter-cluster interference for user ', num2str(i), ': ', num2str(inteer_i(k, i))]);
%     %     disp(['Backscatter inter-cluster interference for user ', num2str(k), ': ', num2str(inteer_b(k, i))]);
%     %     disp(['Total inter-cluster interference for user ', num2str(k), ': ', num2str(inteer_b_all(k, i))]);
%     % end

   


%     % disp(obj_prev)

    

%     % [V_opt, A_opt, B_opt, Ac_opt, Bc_opt, obj_curr, status] = optimize_feasibility_pbf(para,W_init,channel_data,decoding_order,...
%     %  A, B, Ac, Bc,J_r,J_t);

%     %  [V_opt, A_opt, B_opt, A_c_opt, B_c_opt, obj, status] =sca_rate_max_pbf_init(para,W,channel_data,decoding_order,...
%     % A, B, Ac, Bc,channel_data.alpha,J_r,J_t);

%     %  disp(A_opt);
%     %  disp(Ac_opt);
%     %  disp(Bc_opt);
%     %  disp(Bc_opt);

%     %  disp(obj);
%     %  disp(status);


    
% %     %% ITERATIVE OPTIMIZATION
% %     % Stage 1: Feasibility optimization
% %     [W_opt, A, B, Ac, Bc, obj_history] = optimize_feasibility(...
% %         para, H_eff, Hc_eff, channel_data.alpha, A, B, Ac, Bc, decoding_order, W_init);
    
% %     % Extract beamforming vectors
% %     W = extract_beamforming_vectors(W_opt);

% %     % DRIS: Identity J matrices
% %     J_r = eye(N);
% %     J_t = eye(N);

    
% %    [V_opt, A_opt, B_opt, A_c_opt, B_c_opt, obj_prev, status] = feasible_points_passive(para,W,channel_data,decoding_order,...
% %     A, B, Ac, Bc,channel_data.alpha,J_r,J_t);
    
% %     % Stage 2: Sum-rate maximization
% %     [W_opt, A, B, Ac, Bc, rate_history] = maximize_sum_rate_s(...
% %         para, H_eff, Hc_eff, channel_data.alpha, A, B, Ac, Bc, decoding_order, W);
    
% %     % Store results for this Monte Carlo iteration
% %     store_results(mc, W_opt, A, B, Ac, Bc, rate_history);
