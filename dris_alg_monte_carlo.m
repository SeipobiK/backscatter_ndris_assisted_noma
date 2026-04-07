clear all; clc;
addpath(genpath(pwd));
rng(2022); % For reproducibility
diary('results_passive/output_log_test.txt')
diary on
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

alpha_f_mc_dris = zeros(K, MC_MAX);
alpha_n_mc_dris = zeros(K, MC_MAX);

rng_seeds = randi(1e6, MC_MAX, 1);

% Create results directory
results_dir = 'results_passive';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% % % Start parallel pool
if isempty(gcp('nocreate'))
    num_workers = 14;
    pool = parpool('local', num_workers);
    fprintf('Using %d workers for parallel processing\n', pool.NumWorkers);
else
    pool = gcp;
    fprintf('Existing pool with %d workers found\n', pool.NumWorkers);
end

tic;
[BS_array_par, RIS_array_par] = generate_arrays(para);

%% MAIN MONTE CARLO LOOP
fprintf('Starting Monte Carlo iterations...\n');
num_mc_iterations=1;

parfor mc = 1:para.MC_MAX
    try
        fprintf('Monte Carlo Iteration %d/%d\n', mc, num_mc_iterations);
        
        % Set random seed for reproducibility
        rng(rng_seeds(mc), 'twister');
        
        % Generate channels
        [H_local, g_local, f_local] = generate_channel(para, BS_array_par, RIS_array_par);
        
        % Preprocess channels
        channel_data = preprocess_channels(para, H_local, g_local, f_local);
        
        %% Run DRIS (J matrices are identity)
        J_r_dris = eye(N);
        J_t_dris = eye(N);
        
        [Rates_dris, obj_history_dris_mc, w_k_dris_mc, theta_dris_mc, ...
         rate_f_dris, rate_n_dris, rate_c_dris, ...
         noma_signal_dris_mc, BST_signal_dris_mc, ...
         noma_interference_dris_mc, BST_interference_dris_mc, ...
         intra_cluster_dris_mc, inter_cluster_dris_mc, ...
         inter_cluster_BST_dris_mc, inter_cluster_BST_all_dris_mc, ...
         decoding_order_dris_mc, alpha_f_dris, alpha_n_dris] = ...
         run_optimization(para, channel_data, J_r_dris, J_t_dris);
        %  disp(size(obj_history_dris_mc));
        
        % Store DRIS results
        obj_history_dris(:, mc) = obj_history_dris_mc;
        w_k_dris(:, :, :, mc) = w_k_dris_mc;
        theta_dris(:, :, :, mc) = theta_dris_mc;
        rate_f_mc_dris_outer(:, :, mc) = rate_f_dris;
        rate_n_mc_dris_outer(:, :, mc) = rate_n_dris;
        rate_c_mc_dris_outer(:, :, mc) = rate_c_dris;
        rates_dris(:, :, :, mc) = Rates_dris;
        
        noma_signal_dris(:, :, :, mc) = noma_signal_dris_mc;
        BST_signal_dris(:, :, mc) = BST_signal_dris_mc;
        noma_interference_dris(:, :, :, mc) = noma_interference_dris_mc;
        BST_interference_dris(:, :, mc) = BST_interference_dris_mc;
        intra_cluster_dris(:, :, :, mc) = intra_cluster_dris_mc;
        inter_cluster_dris(:, :, :, mc) = inter_cluster_dris_mc;
        inter_cluster_BST_dris(:, :, :, mc) = inter_cluster_BST_dris_mc;
        inter_cluster_BST_all_dris(:, :, :, mc) = inter_cluster_BST_all_dris_mc;
        decoding_order_dris(:, :, :, mc) = decoding_order_dris_mc;
        
        alpha_f_mc_dris(:, mc) = alpha_f_dris;
        alpha_n_mc_dris(:, mc) = alpha_n_dris;
        
        % disp(['MC ', num2str(mc), ' completed - Final DRIS WSR: ', num2str(obj_history_dris_mc(end))]);
        % disp('Rates for DRIS: ');
        % disp(Rates_dris);
        % disp('Objective history for DRIS: ');
        % disp(obj_history_dris_mc);
    catch ME
        fprintf('Error in MC run %d: %s\n', mc, ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s at line %d\n', ME.stack(i).name, ME.stack(i).line);
        end
        obj_history_dris(:, mc) = NaN;
    end
end
% diary off 
toc;

% Calculate averages (excluding NaN runs)
valid_dris = find(~isnan(obj_history_dris(1, :)));
if ~isempty(valid_dris)
    avg_dris = mean(obj_history_dris(:, valid_dris), 2);
    std_dris = std(obj_history_dris(:, valid_dris), 0, 2);
else
    avg_dris = NaN(outer_iter+1, 1);
    std_dris = NaN(outer_iter+1, 1);
    fprintf('Warning: No valid DRIS runs\n');
end

% Plot results with error bars
fig = figure;
x = 0:outer_iter;

% DRIS plot with error bars
plot(x, avg_dris, '-s', 'DisplayName', 'DRIS', ...
    'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.2, 0.2, 0.8]);
hold off;

legend('Location', 'southeast');
xlabel('Iteration Number');
ylabel('Average Weighted Sum Rate (bps/Hz)');
title(sprintf('DRIS Performance: M=%d, N=%d, K=%d, K_c=%d', M, N, K, K_c));
grid on;
xlim([0, outer_iter]);
xticks(0:outer_iter);

% ============================
% Generate timestamp + filenames
% ============================
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

filename = sprintf('dris_results_M%d_N%d_K%d_%s.mat', M, N, K, timestamp);
filename_all = sprintf('full_workspace_dris_M%d_N%d_%s.mat', M, N, timestamp);
filename_plot = sprintf('convergence_dris_M%d_N%d_%s.png', M, N, timestamp);

% Save results structure
results.obj_history_dris = obj_history_dris;
results.avg_dris = avg_dris;
results.std_dris = std_dris;
results.para = para;
results.valid_dris = valid_dris;

% Also save interference metrics if needed
results.rates_dris = rates_dris;
results.noma_signal_dris = noma_signal_dris;
results.BST_signal_dris = BST_signal_dris;
results.noma_interference_dris = noma_interference_dris;
results.BST_interference_dris = BST_interference_dris;
results.intra_cluster_dris = intra_cluster_dris;
results.inter_cluster_dris = inter_cluster_dris;
results.inter_cluster_BST_dris = inter_cluster_BST_dris;
results.inter_cluster_BST_all_dris = inter_cluster_BST_all_dris;
results.decoding_order_dris = decoding_order_dris;
results.alpha_f_mc_dris = alpha_f_mc_dris;
results.alpha_n_mc_dris = alpha_n_mc_dris;
results.w_k_dris = w_k_dris;
results.theta_dris = theta_dris;

% ============================
% Build folder structure
% ============================
base_folder = 'results';
year_str = datestr(now, 'yyyy');
month_str = lower(datestr(now, 'mmm'));
day_str = datestr(now, 'dd');
output_folder = fullfile(base_folder, year_str, month_str, day_str);

% Create folder if it does not exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% ============================
% Save files inside the folder
% ============================
save(fullfile(output_folder, filename_all), '-v7.3');
save(fullfile(output_folder, filename), 'results');
saveas(fig, fullfile(output_folder, filename_plot));

fprintf('\nResults saved in: %s\n', fullfile(output_folder, filename));
fprintf('Workspace saved in: %s\n', fullfile(output_folder, filename_all));
fprintf('Plot saved in: %s\n', fullfile(output_folder, filename_plot));
fprintf('Valid DRIS runs: %d/%d\n', length(valid_dris), MC_MAX);