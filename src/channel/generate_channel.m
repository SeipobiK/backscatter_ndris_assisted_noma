function [H, g, f] = generate_channel(para, BS_array, RIS_array)
    % Rician factor
    epsilon = para.rician;

    BS_loc = para.BS_loc;
    userloc = para.userloc;   % 3 x Nclusters x 3 (d,elev,az)
    Nclusters = para.K;       % number of clusters
    Kusers = para.K_u;        % users per cluster for RIS-user channel (3)
    
    % ============================
    % 1. BS -> RIS channel
    % ============================
    H_NLOS = 1/sqrt(2) * (randn(para.N, para.M) + 1i*randn(para.N, para.M));
    a_BR = steering_vector(BS_array, -BS_loc(2), -BS_loc(3));
    a_RB = steering_vector(RIS_array, BS_loc(2), BS_loc(3));
    H_LOS = a_RB * a_BR.';
   
    path_loss = sqrt(10.^(-para.pathloss(BS_loc(1))/10));
    H = path_loss * ( sqrt(epsilon/(epsilon+1)) * H_LOS + sqrt(1/(epsilon+1)) * H_NLOS );

    % ============================
    % 2. RIS -> Users channel
    % ============================
    g = zeros(para.N, Nclusters, Kusers);
    g_NLOS = 1/sqrt(2) * (randn(para.N, Nclusters, Kusers) + 1i*randn(para.N, Nclusters, Kusers));
    % disp(['Norm of NLOS component for RIS-User channel (before path loss): ', num2str((norm(g_NLOS(:)))^2)]);
    g_LOS = zeros(para.N, Nclusters, Kusers);

    for c = 1:Nclusters
        for k = 1:Kusers
            g_LOS(:,c,k) = steering_vector(RIS_array, userloc(k,c,2), userloc(k,c,3));
            d = userloc(k,c,1); 
            pl = sqrt(10.^(-para.pathloss(d)/10));
            disp(['User ', num2str(k), ' in Cluster ', num2str(c), ' distance: ', num2str(d), ' m, path loss: ', num2str(pl), ...
            ' Norm of LOS component: ', num2str(norm(g_LOS(:,c,k))^2/para.N),' Norm of NLOS component: ', num2str(norm(g_NLOS(:,c,k))^2/para.N)]);
            g(:,c,k) = pl * ( sqrt(epsilon/(epsilon+1)) * g_LOS(:,c,k) + sqrt(1/(epsilon+1)) * g_NLOS(:,c,k) );
        end
    end

    % ============================
    % 3. BS -> Users channel (f)
    % Generalized for any number of users per cluster (K_c)
    % ============================
    
    % Determine number of users per cluster for BS-User channel
    % If para.K_c is defined, use it; otherwise use Kusers (for backward compatibility)
    if isfield(para, 'K_c')
        K_c = para.K_c;
    else
        K_c = Kusers; % Default to RIS-user channel user count
    end
    
    f = zeros(Nclusters, K_c); % BS->User channels: Nclusters x K_c
    
    for c = 1:Nclusters
        % Get distances for all users in this cluster
        cluster_dists = zeros(1, K_c);
        
        % Check if distances are stored in para.BSTdist structure
        if isfield(para, 'BSTdist') && isfield(para.BSTdist, ['Cluster' num2str(c)])
            % Method 1: Using stored distances from para.BSTdist
            for k = 1:K_c
                user_field = ['User' num2str(k) '_BD'];
                if isfield(para.BSTdist.(['Cluster' num2str(c)]), user_field)
                    cluster_dists(k) = para.BSTdist.(['Cluster' num2str(c)]).(user_field);
                else
                    % Fallback: use distances from userloc if available
                    if k <= size(userloc, 1)
                        cluster_dists(k) = userloc(k, c, 1);
                    else
                        % Default distance if not specified
                        cluster_dists(k) = para.default_BS_user_distance;
                    end
                end
            end
        elseif size(userloc, 1) >= K_c
            % Method 2: Use distances from userloc array
            for k = 1:K_c
                cluster_dists(k) = userloc(k, c, 1);
            end
        else
            % Method 3: Generate random distances within a reasonable range
            for k = 1:K_c
                % Random distance between 20m and 200m
                cluster_dists(k) = 20 + 180 * rand();
            end
        end
        
        % Generate channel for each user in the cluster
        for k = 1:K_c
            dist = cluster_dists(k);
            
            % Path loss (convert from dB to linear)
            pl = sqrt(10^(-para.pathloss(dist)/10));
            
            % Small-scale Rayleigh fading (MIMO channel: M antennas)
            % Use MIMO Rayleigh fading (M x 1 vector for each user)
            f(c, k) = pl * (randn + 1i*randn) / sqrt(2);
        end
    end
end