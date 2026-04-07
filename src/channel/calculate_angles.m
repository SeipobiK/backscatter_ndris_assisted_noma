function [angles, users,distances] = calculate_angles()
    % ============================
    % Fixed BS and RIS positions
    % ============================
BS = [0, 10, 20];
RIS = [0, 30, 20];

% ============================
% Hard-coded cluster centers
% ============================
cluster_centers = [0, 25, 0;
                   0, 35, 0;
                   5, 35, 0]; % new cluster center

radii = [5,5,5]; % radius for each cluster

% ============================
% Hard-coded users per cluster (3 users each)
% ============================

%     % ----- Cluster 1 -----
    % BD is at index 5 (last user)
    users{1} = [ 2  , 28, 0;   % NOMA user 1
                 1, 26, 0;  % NOMA user 2
                 -2, 25, 0;  % NOMA user 3
                 2, 24, 0;   % NOMA user 4
                 2, 27, 0];  % BD (Backscatter Device)

    % ----- Cluster 2 -----
    % BD is at index 5 (last user)
    users{2} = [ 2, 30, 0;   % NOMA user 1
                 4, 30, 0;   % NOMA user 2
                 6, 30, 0;   % NOMA user 3
                 7, 30, 0;   % NOMA user 4
                 4, 29, 0];  % BD (Backscatter Device)



% % ----- Cluster 1 (center: [0,25,0], radius: 5m) -----
% users{1} = [ 2, 28, 0;   
%              1, 26, 0;   
%              2, 27, 0];  

% % ----- Cluster 2 (center: [0,35,0], radius: 5m) -----
% % users{2} = [ -1, 31, 0;   % User 1: North-East
% %              -1, 40, 0;   % User 2: South-West
% %              0, 31, 0];   % User 3: North-East
% % 
% % ----- Cluster 3 (center: [5,30,0], radius: 5m) -----
% users{2} = [ 4, 31, 0;   % User 1: North-East   
%              6, 30, 0;   % User 2: South-West
%              4, 29, 0];  % User 3: South-East



% ----- Define cluster centers and radii -----
cluster_centers = [0 25 0;   % Cluster 1 center
                   0 35 0;   % Cluster 2 center
                   5 30 0];  % Cluster 3 center
    radius = 5;                  % Cluster radius (meters)
    
    num_clusters = size(cluster_centers,1);
    num_users_per_cluster = 3;
    
    % ----- Generate user positions -----
    % users = cell(num_clusters,1);
    % ----- Cluster 1 (center: [0,25,0], away direction +Y) -----
% users{1} = [ -1, 20, 0;     % left
%               0,  20, 0;     % center
%               1,  20, 0];    % right
% 
% % ----- Cluster 2 (center: [0,35,0], away direction +Y) -----
% users{2} = [ -1, 40, 0;
%               0,  40, 0;
%               1,  40, 0];
% 
% % ----- Cluster 3 (center: [5,30,0], away direction +X+Y) -----
% users{3} = [ 10, 31, 0;
%               9, 32, 0;
%               7, 26, 0];

    Nclusters = length(users);
    disp(Nclusters);
    Kusers = size(users{1},1);

    % ============================
    % 1. AoD at BS (BS -> RIS)
    % ============================
    delta_BS_RIS = RIS - BS;
    angles.BS_AoD.azimuth = atan2d(delta_BS_RIS(2), delta_BS_RIS(1));
    angles.BS_AoD.elevation = atan2d(delta_BS_RIS(3), sqrt(delta_BS_RIS(1)^2 + delta_BS_RIS(2)^2));
    angles.BS_AoD.distance = norm(delta_BS_RIS);

    % ============================
    % 2. AoA at RIS
    % ============================
    angles.RIS_AoA.azimuth = mod(angles.BS_AoD.azimuth + 180, 360);
    angles.RIS_AoA.elevation = -angles.BS_AoD.elevation;
    angles.RIS_AoA.distance = angles.BS_AoD.distance;

    % ============================
    % 3. AoD at RIS -> Users
    % ============================
    for c = 1:Nclusters
        for k = 1:Kusers
            user = users{c}(k,:);
            delta_RIS_user = user - RIS;

            azimuth = atan2d(delta_RIS_user(2), delta_RIS_user(1));
            elevation = atan2d(delta_RIS_user(3), sqrt(delta_RIS_user(1)^2 + delta_RIS_user(2)^2));

            angles.RIS_AoD.clusters(c).users(k).user = sprintf('C%d-U%d', c, k);
            angles.RIS_AoD.clusters(c).users(k).azimuth = azimuth;
            angles.RIS_AoD.clusters(c).users(k).elevation = elevation;
            angles.RIS_AoD.clusters(c).users(k).distance = norm(delta_RIS_user);
        end
    end
    
    % Initialize struct array
    distances = struct();
   
    
    for c = 1:length(users)
        BD = users{c}(3, :); % BD position
        clusterDistances = struct();
        
        % Distance between each user and BD
        for u = 1:2
            user_pos = users{c}(u, :);
            d_user_BD = norm(user_pos - BD);
            clusterDistances.(['User' num2str(u) '_BD']) = d_user_BD;
        end
        
        % Optional: distance between users
        d_user1_user2 = norm(users{c}(1,:) - users{c}(2,:));
        clusterDistances.User1_User2 = d_user1_user2;
    
        % Store in main struct
        distances.(['Cluster' num2str(c)]) = clusterDistances;
    end

    % 
%     % % % % % % % % % ============================
%     % % % % % % % % 4. Visualization
%     % % % % % % % ============================
    figure; hold on; grid on; axis equal;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('BS, RIS, and User Clusters with Circles');

    % Plot BS and RIS
    scatter3(BS(1),BS(2),BS(3),100,'r','filled','^'); text(BS(1),BS(2),BS(3)+1,'BS');
    scatter3(RIS(1),RIS(2),RIS(3),100,'b','filled','s'); text(RIS(1),RIS(2),RIS(3)+1,'RIS');

    colors = lines(Nclusters);
    for c = 1:Nclusters
        % Cluster center
        center = cluster_centers(c,:);
        scatter3(center(1), center(2), center(3), 80, colors(c,:),'d','filled');
        text(center(1), center(2), center(3)+0.5,sprintf('Cluster %d', c),'Color','k','FontSize',10,'FontWeight','bold');

        % Draw cluster circle in XY-plane
        theta = linspace(0, 2*pi, 100);
        x_circle = center(1) + radii(c) * cos(theta);
        y_circle = center(2) + radii(c) * sin(theta);
        z_circle = center(3) * ones(size(theta));
        plot3(x_circle, y_circle, z_circle, '--','Color', colors(c,:), 'LineWidth',1.2);

        % Users
        cluster_users = users{c};
        scatter3(cluster_users(:,1), cluster_users(:,2), cluster_users(:,3), 50, colors(c,:),'filled','o');

        % Label users
        for k = 1:Kusers
            text(cluster_users(k,1), cluster_users(k,2), cluster_users(k,3)+0.5, ...
                sprintf('C%d-U%d', c, k), 'Color', colors(c,:), 'FontSize', 8, 'FontWeight','bold');
        end
    end
    legend({'BS','RIS','Cluster Centers','Cluster Circles','Users'});
    view(45,25);
    hold off;
end



% function [angles, users, distances] = calculate_angles()
%     % ============================
%     % Fixed BS and RIS positions
%     % ============================
%     BS = [0, 10, 20];
%     RIS = [0, 30, 20];

%     % ============================
%     % Hard-coded users per cluster
%     % Last user in each cluster is the BD
%     % First N-1 users are NOMA users
%     % ============================

%     % ----- Cluster 1 -----
%     % BD is at index 5 (last user)
%     users{1} = [ 2  , 28, 0;   % NOMA user 1
%                  1, 26, 0;  % NOMA user 2
%                  -2, 25, 0;  % NOMA user 3
%                  2, 24, 0;   % NOMA user 4
%                  2, 27, 0];  % BD (Backscatter Device)

%     % ----- Cluster 2 -----
%     % BD is at index 5 (last user)
%     users{2} = [ 2, 31, 0;   % NOMA user 1
%                  6, 30, 0;   % NOMA user 2
%                  6, 30, 0;   % NOMA user 3
%                  7, 30, 0;   % NOMA user 4
%                  4, 29, 0];  % BD (Backscatter Device)

%     % ============================
%     % Cluster parameters for visualization
%     % ============================
%     cluster_centers = [0, 25, 0;    % Cluster 1 center
%                        5, 30, 0];   % Cluster 2 center
%     radii = [5, 5];                  % Radius for each cluster

%     Nclusters = length(users);
    
%     % Get number of NOMA users per cluster (all except last/BD)
%     for c = 1:Nclusters
%         num_users_per_cluster(c) = size(users{c}, 1);
%         num_noma_users(c) = num_users_per_cluster(c) - 1;  % Last user is BD
%     end

%     % ============================
%     % 1. AoD at BS (BS -> RIS)
%     % ============================
%     delta_BS_RIS = RIS - BS;
%     angles.BS_AoD.azimuth = atan2d(delta_BS_RIS(2), delta_BS_RIS(1));
%     angles.BS_AoD.elevation = atan2d(delta_BS_RIS(3), sqrt(delta_BS_RIS(1)^2 + delta_BS_RIS(2)^2));
%     angles.BS_AoD.distance = norm(delta_BS_RIS);

%     % ============================
%     % 2. AoA at RIS (from BS)
%     % ============================
%     angles.RIS_AoA.azimuth = mod(angles.BS_AoD.azimuth + 180, 360);
%     angles.RIS_AoA.elevation = -angles.BS_AoD.elevation;
%     angles.RIS_AoA.distance = angles.BS_AoD.distance;

%     % ============================
%     % 3. AoD at RIS -> Users (NOMA users + BD)
%     % ============================
%     for c = 1:Nclusters
%         num_users = size(users{c}, 1);
        
%         for k = 1:num_users
%             user = users{c}(k, :);
%             delta_RIS_user = user - RIS;

%             azimuth = atan2d(delta_RIS_user(2), delta_RIS_user(1));
%             elevation = atan2d(delta_RIS_user(3), sqrt(delta_RIS_user(1)^2 + delta_RIS_user(2)^2));

%             % Determine if user is NOMA or BD
%             if k < num_users
%                 user_type = 'NOMA';
%             else
%                 user_type = 'BD';
%             end
            
%             angles.RIS_AoD.clusters(c).users(k).user_id = sprintf('C%d-%s%d', c, user_type, k);
%             angles.RIS_AoD.clusters(c).users(k).type = user_type;
%             angles.RIS_AoD.clusters(c).users(k).azimuth = azimuth;
%             angles.RIS_AoD.clusters(c).users(k).elevation = elevation;
%             angles.RIS_AoD.clusters(c).users(k).distance = norm(delta_RIS_user);
%             angles.RIS_AoD.clusters(c).users(k).position = user;
%         end
%     end

%     % ============================
%     % 4. Distance Calculations
%     % ============================
%     distances = struct();
    
%     for c = 1:Nclusters
%         cluster_users = users{c};
%         num_users = size(cluster_users, 1);
        
%         % Last user is BD
%         BD_position = cluster_users(end, :);
%         noma_positions = cluster_users(1:end-1, :);
%         num_noma = size(noma_positions, 1);
        
%         clusterDistances = struct();
        
%         % Distance between each NOMA user and BD
%         for u = 1:num_noma
%             d_noma_bd = norm(noma_positions(u, :) - BD_position);
%             clusterDistances.(sprintf('NOMA%d_to_BD', u)) = d_noma_bd;
%         end
        
%         % Distances between NOMA users (for NOMA pairing)
%         noma_distances = zeros(num_noma);
%         for u1 = 1:num_noma
%             for u2 = u1+1:num_noma
%                 d_noma_noma = norm(noma_positions(u1, :) - noma_positions(u2, :));
%                 clusterDistances.(sprintf('NOMA%d_NOMA%d', u1, u2)) = d_noma_noma;
%                 noma_distances(u1, u2) = d_noma_noma;
%             end
%         end
%         clusterDistances.noma_pairwise_matrix = noma_distances;
        
%         % Distance from RIS to BD and NOMA users
%         clusterDistances.RIS_to_BD = norm(BD_position - RIS);
%         clusterDistances.RIS_to_NOMA = zeros(1, num_noma);
%         for u = 1:num_noma
%             clusterDistances.RIS_to_NOMA(u) = norm(noma_positions(u, :) - RIS);
%         end
        
%         % Distance from BS to BD and NOMA users
%         clusterDistances.BS_to_BD = norm(BD_position - BS);
%         clusterDistances.BS_to_NOMA = zeros(1, num_noma);
%         for u = 1:num_noma
%             clusterDistances.BS_to_NOMA(u) = norm(noma_positions(u, :) - BS);
%         end
        
%         % Store in main struct
%         distances.(['Cluster' num2str(c)]) = clusterDistances;
%     end

%     % ============================
%     % 5. Visualization
%     % ============================
%     figure; hold on; grid on; axis equal;
%     xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
%     title('BS, RIS, and User Clusters with NOMA Users and BD');

%     % Plot BS and RIS
%     scatter3(BS(1), BS(2), BS(3), 100, 'r', 'filled', '^'); 
%     text(BS(1), BS(2), BS(3)+1, 'BS', 'FontSize', 10, 'FontWeight', 'bold');
    
%     scatter3(RIS(1), RIS(2), RIS(3), 100, 'b', 'filled', 's'); 
%     text(RIS(1), RIS(2), RIS(3)+1, 'RIS', 'FontSize', 10, 'FontWeight', 'bold');

%     colors = lines(Nclusters);
    
%     for c = 1:Nclusters
%         % Cluster center
%         center = cluster_centers(c, :);
%         scatter3(center(1), center(2), center(3), 80, colors(c, :), 'd', 'filled');
%         text(center(1), center(2), center(3)+0.5, sprintf('Cluster %d', c), ...
%             'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');

%         % Draw cluster circle in XY-plane
%         theta = linspace(0, 2*pi, 100);
%         x_circle = center(1) + radii(c) * cos(theta);
%         y_circle = center(2) + radii(c) * sin(theta);
%         z_circle = center(3) * ones(size(theta));
%         plot3(x_circle, y_circle, z_circle, '--', 'Color', colors(c, :), 'LineWidth', 1.2);

%         % Users
%         cluster_users = users{c};
%         num_users = size(cluster_users, 1);
        
%         % Plot NOMA users (all except last)
%         noma_users = cluster_users(1:end-1, :);
%         scatter3(noma_users(:,1), noma_users(:,2), noma_users(:,3), 50, colors(c, :), 'filled', 'o');
        
%         % Plot BD (last user) with different marker
%         BD = cluster_users(end, :);
%         scatter3(BD(1), BD(2), BD(3), 80, colors(c, :), 'filled', 'diamond');
        
%         % Label users
%         for k = 1:num_users-1  % NOMA users
%             text(cluster_users(k,1), cluster_users(k,2), cluster_users(k,3)+0.5, ...
%                 sprintf('NOMA%d', k), 'Color', colors(c, :), 'FontSize', 8, 'FontWeight', 'bold');
%         end
%         % Label BD
%         text(BD(1), BD(2), BD(3)+0.5, sprintf('BD'), ...
%             'Color', colors(c, :), 'FontSize', 9, 'FontWeight', 'bold');
%     end
    
%     legend({'BS', 'RIS', 'Cluster Centers', 'Cluster Circles', 'NOMA Users', 'BD'}, ...
%         'Location', 'best');
%     view(45, 25);
%     hold off;
    
%     % Display summary
%     fprintf('\n=== Network Configuration Summary ===\n');
%     for c = 1:Nclusters
%         fprintf('Cluster %d: %d NOMA users, 1 BD\n', c, num_noma_users(c));
%         fprintf('  BS to RIS distance: %.2f m\n', angles.BS_AoD.distance);
%         fprintf('  RIS to BD distance: %.2f m\n', distances.(['Cluster' num2str(c)]).RIS_to_BD);
%         fprintf('  BS to BD distance: %.2f m\n', distances.(['Cluster' num2str(c)]).BS_to_BD);
%     end
% end