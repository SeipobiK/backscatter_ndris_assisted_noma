clear all; clc;
addpath(genpath(pwd));
rng(2022); % For reproducibility
% diary('results_passive/output_log_test.txt')
% diary on
% Initialize parameters
para = para_init();
[BS_array, RIS_array] = generate_arrays(para);

% Constants
outer_iter = para.outer_iter;
MC_MAX = para.MC_MAX;
K = para.K;
K_c = para.K_c;
N = para.N;
M = para.M;
num_mc_iterations = MC_MAX;

% Preallocate results for DRIS
obj_history_dris = zeros(outer_iter+1, MC_MAX);
w_k_dris = zeros(M, K, outer_iter+1, MC_MAX);
theta_dris = zeros(N, N, outer_iter+1, MC_MAX);

rate_f_mc_dris_outer = zeros(K, outer_iter+1, MC_MAX);
rate_n_mc_dris_outer = zeros(K, outer_iter+1, MC_MAX);
rate_c_mc_dris_outer = zeros(K, outer_iter+1, MC_MAX);
rates_dris = zeros(K, K_c, outer_iter+1, MC_MAX);

% Interference and signal metrics for DRIS
noma_signal_dris = zeros(K, K_c, outer_iter+1, MC_MAX);
BST_signal_dris = zeros(K, outer_iter+1, MC_MAX);
noma_interference_dris = zeros(K, K_c, outer_iter+1, MC_MAX);
BST_interference_dris = zeros(K, outer_iter+1, MC_MAX);
intra_cluster_dris = zeros(K, K_c, outer_iter+1, MC_MAX);
inter_cluster_dris = zeros(K, K_c, outer_iter+1, MC_MAX);
inter_cluster_BST_dris = zeros(K, K_c, outer_iter+1, MC_MAX);
inter_cluster_BST_all_dris = zeros(K, K_c, outer_iter+1, MC_MAX);
decoding_order_dris = zeros(K, K_c, outer_iter+1, MC_MAX);

alpha_history_dris = zeros(K, K_c, outer_iter+1, MC_MAX);
gains_it_history_dris = zeros(K, K_c, outer_iter+1, MC_MAX);


rng_seeds = randi(1e6, MC_MAX, 1);

% % Create results directory
results_dir = 'results_passive';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% % % Start parallel pool
% if isempty(gcp('nocreate'))
%     num_workers = 14;
%     pool = parpool('local', num_workers);
%     fprintf('Using %d workers for parallel processing\n', pool.NumWorkers);
% else
%     pool = gcp;
%     fprintf('Existing pool with %d workers found\n', pool.NumWorkers);
% end

tic;
[BS_array_par, RIS_array_par] = generate_arrays(para);

%% MAIN MONTE CARLO LOOP
fprintf('Starting Monte Carlo iterations...\n');
num_mc_iterations=1;

for mc = 10:para.MC_MAX
    % try
        fprintf('Monte Carlo Iteration %d/%d\n', mc, num_mc_iterations);
        
        % Set random seed for reproducibility
        rng(rng_seeds(mc), 'twister');

        % disp(alpha);
        
        % Generate channels
        [H_local, g_local, f_local] = generate_channel(para, BS_array_par, RIS_array_par);
        
        % Preprocess channels
        channel_data = preprocess_channels(para, H_local, g_local, f_local);
        
        %% Run DRIS (J matrices are identity)
        J_r_dris = eye(N);
        J_t_dris = eye(N);
        
        [Rates,obj_history, w_k_history, theta_history, ...
          rate_f_history, rate_n_history, rate_c_history, ...
          noma_signal_history, BST_signal_history, ...
          noma_interference_history, BST_interference_history, ...
          intra_cluster_history, inter_cluster_history, ...
          inter_cluster_BST_history, inter_cluster_BST_all_history, ...
          decoding_order_history, alpha_f, alpha_n] = channel_verification(para, channel_data, J_r_dris, J_t_dris);
        %  disp(size(obj_history_dris_mc));
        
        % Store DRIS results
        obj_history_dris(:, mc) = obj_history;
        % w_k_dris(:, :, :, mc) = w_k_dris_mc;
        % theta_dris(:, :, :, mc) = theta_dris_mc;
        % rate_f_mc_dris_outer(:, :, mc) = rate_f_dris;
        % rate_n_mc_dris_outer(:, :, mc) = rate_n_dris;
        % rate_c_mc_dris_outer(:, :, mc) = rate_c_dris;
        % rates_dris(:, :, :, mc) = Rates_dris;
        
        % noma_signal_dris(:, :, :, mc) = noma_signal_dris_mc;
        % BST_signal_dris(:, :, mc) = BST_signal_dris_mc;
        % noma_interference_dris(:, :, :, mc) = noma_interference_dris_mc;
        % BST_interference_dris(:, :, mc) = BST_interference_dris_mc;
        % intra_cluster_dris(:, :, :, mc) = intra_cluster_dris_mc;
        % inter_cluster_dris(:, :, :, mc) = inter_cluster_dris_mc;
        % inter_cluster_BST_dris(:, :, :, mc) = inter_cluster_BST_dris_mc;
        % inter_cluster_BST_all_dris(:, :, :, mc) = inter_cluster_BST_all_dris_mc;
        % decoding_order_dris(:, :, :, mc) = decoding_order_dris_mc;
        % alpha_history_dris(:, :, :, mc) = alpha_history_dris_mc;
        % gains_it_history_dris(:, :, :, mc) = gains_it_history_dris_mc;

        
        % alpha_f_mc_dris(:, mc) = alpha_f_dris;
        % alpha_n_mc_dris(:, mc) = alpha_n_dris;
        
        % disp(['MC ', num2str(mc), ' completed - Final DRIS WSR: ', num2str(obj_history_dris_mc(end))]);
        % disp('Rates for DRIS: ');
        % disp(Rates_dris);
        % disp('Objective history for DRIS: ');
        % disp(obj_history_dris_mc);
    % catch ME
    %     fprintf('Error in MC run %d: %s\n', mc, ME.message);
    %     fprintf('Stack trace:\n');
    %     for i = 1:length(ME.stack)
    %         fprintf('  %s at line %d\n', ME.stack(i).name, ME.stack(i).line);
    %     end
    %     obj_history_dris(:, mc) = NaN;
    % end
end