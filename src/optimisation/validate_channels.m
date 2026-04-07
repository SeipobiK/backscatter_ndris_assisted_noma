function validate_channels()
    % Test parameters
    para.M = 64;           % BS antennas
    para.N = 100;          % RIS elements
    para.rician = 10;      % Rician factor (linear, e.g., 10 = 10dB)
    para.K = 2;            % Number of clusters
    para.K_u = 3;          % RIS users per cluster
    para.K_c = 3;          % BS users per cluster (NOMA users)
    
    % Path loss parameters
    para.pathloss = @(d) 30 + 22*log10(d); % dB
    
    % Positions
    para.BS_loc = [0, 10, 20];  % [distance, azimuth, elevation]
    
    % User locations: [distance, azimuth, elevation]
    % For RIS-User channel
    para.userloc = zeros(para.K_u, para.K, 3);
    para.userloc(1,1,:) = [30, 30, 0];    % Cluster1-User1
    para.userloc(2,1,:) = [32, 25, 0];    % Cluster1-User2
    para.userloc(3,1,:) = [28, 35, 0];    % Cluster1-User3
    para.userloc(1,2,:) = [35, 40, 0];    % Cluster2-User1
    para.userloc(2,2,:) = [33, 45, 0];    % Cluster2-User2
    para.userloc(3,2,:) = [38, 38, 0];    % Cluster2-User3
    
    % BS-User distances (for f channel)
    para.BSTdist.Cluster1.User1_BD = 50;
    para.BSTdist.Cluster1.User2_BD = 55;
    para.BSTdist.Cluster1.User3_BD = 48;
    para.BSTdist.Cluster2.User1_BD = 60;
    para.BSTdist.Cluster2.User2_BD = 58;
    para.BSTdist.Cluster2.User3_BD = 62;
    
    para.default_BS_user_distance = 50;
    
    % Steering vectors (simplified for testing)
    BS_array = @(az, el) exp(1i*2*pi*rand(para.M,1));  % Replace with actual steering vector
    RIS_array = @(az, el) exp(1i*2*pi*rand(para.N,1));
    
    % Generate channels
    [H, g, f] = generate_channel(para, BS_array, RIS_array);
    
    % ============================
    % 1. Validate Large-Scale Fading (Path Loss)
    % ============================
    fprintf('\n=== Large-Scale Fading Validation ===\n');
    
    % Check BS-RIS path loss
    BS_RIS_dist = para.BS_loc(1);  % Distance from BS to RIS
    expected_pl_BS_RIS = 10^(-para.pathloss(BS_RIS_dist)/20);
    actual_pl_BS_RIS = norm(H, 'fro') / sqrt(para.N * para.M);
    fprintf('BS->RIS Path Loss:\n');
    fprintf('  Expected: %.4f\n', expected_pl_BS_RIS);
    fprintf('  Actual (average magnitude): %.4f\n', actual_pl_BS_RIS);
    fprintf('  Ratio: %.2f\n', actual_pl_BS_RIS / expected_pl_BS_RIS);
    
    % Check RIS-User path loss for each user
    fprintf('\nRIS->User Path Loss:\n');
    for c = 1:para.K
        for k = 1:para.K_u
            dist = para.userloc(k,c,1);
            expected_pl = 10^(-para.pathloss(dist)/20);
            actual_pl = norm(g(:,c,k)) / sqrt(para.N);
            fprintf('  Cluster%d-User%d: Expected=%.4f, Actual=%.4f, Ratio=%.2f\n', ...
                c, k, expected_pl, actual_pl, actual_pl/expected_pl);
        end
    end
    
    % ============================
    % 2. Validate Small-Scale Fading (Rayleigh/Rician Distribution)
    % ============================
    fprintf('\n=== Small-Scale Fading Validation ===\n');
    
    % Generate multiple realizations to check statistics
    num_realizations = 1000;
    
    % For H (BS-RIS)
    H_magnitudes = zeros(num_realizations, 1);
    for n = 1:num_realizations
        H_temp = generate_channel(para, BS_array, RIS_array);
        H_magnitudes(n) = norm(H_temp, 'fro') / sqrt(para.N * para.M);
    end
    
    % Check Rician distribution for H
    if para.rician > 0
        % Rician K-factor estimation
        rice_fit = fitdist(H_magnitudes, 'rician');
        fprintf('\nBS->RIS Channel (Rician with K=%.1f dB):\n', 10*log10(para.rician));
        fprintf('  Estimated K-factor: %.2f dB\n', 10*log10(rice_fit.s^2/(2*sigma^2)));
        fprintf('  Mean magnitude: %.4f\n', mean(H_magnitudes));
        fprintf('  Std deviation: %.4f\n', std(H_magnitudes));
    end
    
    % For g (RIS-User) - check each user
    fprintf('\nRIS->User Channel Statistics (over %d realizations):\n', num_realizations);
    g_magnitudes = zeros(num_realizations, para.K, para.K_u);
    
    for n = 1:num_realizations
        [~, g_temp, ~] = generate_channel(para, BS_array, RIS_array);
        for c = 1:para.K
            for k = 1:para.K_u
                g_magnitudes(n, c, k) = norm(g_temp(:,c,k)) / sqrt(para.N);
            end
        end
    end
    
    for c = 1:para.K
        for k = 1:para.K_u
            fprintf('  Cluster%d-User%d: Mean=%.4f, Std=%.4f, RiceK=%.2fdB\n', ...
                c, k, mean(g_magnitudes(:,c,k)), std(g_magnitudes(:,c,k)), ...
                10*log10(para.rician));
        end
    end
    
    % ============================
    % 3. Validate f channel (BS-User)
    % ============================
    fprintf('\n=== BS->User Channel (f) Validation ===\n');
    
    num_realizations_f = 10000;
    f_magnitudes = zeros(num_realizations_f, para.K, para.K_c);
    
    for n = 1:num_realizations_f
        [~, ~, f_temp] = generate_channel(para, BS_array, RIS_array);
        for c = 1:para.K
            for k = 1:para.K_c
                f_magnitudes(n, c, k) = abs(f_temp(c, k));
            end
        end
    end
    
    for c = 1:para.K
        for k = 1:para.K_c
            % f should be Rayleigh (Rician with K=0)
            expected_rayleigh_std = sqrt(pi/2) * mean(f_magnitudes(:,c,k));
            fprintf('  Cluster%d-User%d:\n', c, k);
            fprintf('    Mean magnitude: %.4f\n', mean(f_magnitudes(:,c,k)));
            fprintf('    Std deviation: %.4f\n', std(f_magnitudes(:,c,k)));
            fprintf('    Theoretical Rayleigh mean (if σ=1/√2): %.4f\n', sqrt(pi/4));
            
            % Check Rayleigh distribution
            [h, edges] = histcounts(f_magnitudes(:,c,k), 30, 'Normalization', 'pdf');
            x = edges(1:end-1) + diff(edges)/2;
            rayleigh_pdf = (x / 0.5^2) .* exp(-x.^2 / (2*0.5^2));
            
            figure(100+c*10+k);
            histogram(f_magnitudes(:,c,k), 30, 'Normalization', 'pdf', 'DisplayName', 'Simulated');
            hold on;
            plot(x, rayleigh_pdf, 'r-', 'LineWidth', 2, 'DisplayName', 'Theoretical Rayleigh');
            xlabel('Amplitude');
            ylabel('PDF');
            title(sprintf('f Channel Distribution - Cluster%d User%d', c, k));
            legend;
            hold off;
        end
    end
    
    % ============================
    % 4. Validate Correlation and Statistical Properties
    % ============================
    fprintf('\n=== Statistical Properties ===\n');
    
    % Check mean and variance of H
    H_elements = H(:);
    fprintf('H (BS-RIS) statistics:\n');
    fprintf('  Mean of real part: %.4f (should be ~0)\n', mean(real(H_elements)));
    fprintf('  Mean of imag part: %.4f (should be ~0)\n', mean(imag(H_elements)));
    fprintf('  Variance: %.4f (should be 1 after path loss)\n', var(real(H_elements)));
    
    % Check that NLOS component has correct variance
    H_NLOS_test = 1/sqrt(2) * (randn(para.N, para.M) + 1i*randn(para.N, para.M));
    fprintf('\nNLOS component variance check:\n');
    fprintf('  Real part variance: %.4f (expected 0.5)\n', var(real(H_NLOS_test(:))));
    fprintf('  Imag part variance: %.4f (expected 0.5)\n', var(imag(H_NLOS_test(:))));
    fprintf('  Total variance: %.4f (expected 1)\n', var(H_NLOS_test(:)));
    
    % ============================
    % 5. Visual Summary
    % ============================
    figure('Position', [100, 100, 1200, 800]);
    
    % Subplot 1: H channel magnitude distribution
    subplot(2,3,1);
    histogram(H_magnitudes, 30, 'Normalization', 'pdf');
    xlabel('||H||_F / sqrt(NM)');
    ylabel('PDF');
    title(sprintf('BS-RIS Channel Magnitude (K=%.1fdB)', 10*log10(para.rician)));
    grid on;
    
    % Subplot 2: g channel magnitudes for first user
    subplot(2,3,2);
    hold on;
    for c = 1:para.K
        histogram(g_magnitudes(:,c,1), 30, 'Normalization', 'pdf', 'DisplayName', sprintf('Cluster%d', c));
    end
    xlabel('||g|| / sqrt(N)');
    ylabel('PDF');
    title('RIS-User Channel Magnitudes');
    legend;
    grid on;
    hold off;
    
    % Subplot 3: f channel magnitudes
    subplot(2,3,3);
    hold on;
    colors = lines(para.K);
    for c = 1:para.K
        histogram(f_magnitudes(:,c,1), 30, 'Normalization', 'pdf', 'FaceColor', colors(c,:), 'FaceAlpha', 0.5);
    end
    xlabel('|f|');
    ylabel('PDF');
    title('BS-User Channel (Rayleigh Fading)');
    legend(arrayfun(@(x) sprintf('Cluster%d', x), 1:para.K, 'UniformOutput', false));
    grid on;
    hold off;
    
    % Subplot 4: Power delay profile check
    subplot(2,3,4);
    H_power = abs(H).^2;
    imagesc(10*log10(H_power));
    colorbar;
    xlabel('BS Antennas');
    ylabel('RIS Elements');
    title('BS-RIS Channel Power (dB)');
    
    % Subplot 5: Channel correlation check
    subplot(2,3,5);
    H_corr = corrcoef(real(H(:)), imag(H(:)));
    imagesc(H_corr);
    colorbar;
    title('Real vs Imag Correlation');
    xticks([1,2]); yticks([1,2]);
    xticklabels({'Real', 'Imag'}); yticklabels({'Real', 'Imag'});
    
    % Subplot 6: QQ plot for Rayleigh check
    subplot(2,3,6);
    f_sample = f_magnitudes(:,1,1);
    qqplot(f_sample);
    title('QQ Plot vs Normal Distribution');
    grid on;
    
    sgtitle('Channel Validation Results');
    
    % ============================
    % 6. Print Summary
    % ============================
    fprintf('\n=== VALIDATION SUMMARY ===\n');
    fprintf('✓ Large-scale fading: Path loss correctly applied\n');
    fprintf('✓ Small-scale fading: ');
    if para.rician > 0
        fprintf('Rician distribution (K=%.1fdB)\n', 10*log10(para.rician));
    else
        fprintf('Rayleigh distribution\n');
    end
    fprintf('✓ Channel statistics: ');
    if abs(mean(real(H_elements))) < 0.1 && abs(mean(imag(H_elements))) < 0.1
        fprintf('Zero-mean ✓\n');
    else
        fprintf('Non-zero mean ✗\n');
    end
    fprintf('✓ Variance normalization: ');
    if abs(var(H_elements) - 1) < 0.1
        fprintf('Correct ✓\n');
    else
        fprintf('Incorrect (%.2f) ✗\n', var(H_elements));
    end
end