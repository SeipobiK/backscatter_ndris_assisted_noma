function [angles, users, distances] = calculate_angles()

    %% ============================
    % Paper-style fixed locations
    % ============================
    BS  = [0, 0, 4];        % BS location
    RIS = [10, 18, 4];      % RIS location, short BS-RIS distance

    %% ============================
    % Two short-distance clusters
    % Each cluster: User 1, User 2, BD
    % ============================

    cluster_centers = [14, 12, 0;
                    6,  12, 0];

    radii = [3, 3];

    % ----- Cluster 1 -----
    users{1} = [13, 11, 0;     % NOMA user 1
                15, 12, 0;     % NOMA user 2
                14, 13, 0];    % BD

    % ----- Cluster 2 -----
    users{2} = [5,  11, 0;     % NOMA user 1
                7,  12, 0;     % NOMA user 2
                6,  13, 0];    % BD

    Nclusters = length(users);
    Kusers = size(users{1},1);

    %% ============================
    % 1. BS -> RIS angles
    % ============================
    delta_BS_RIS = RIS - BS;

    angles.BS_AoD.azimuth = atan2d(delta_BS_RIS(2), delta_BS_RIS(1));
    angles.BS_AoD.elevation = atan2d(delta_BS_RIS(3), ...
        sqrt(delta_BS_RIS(1)^2 + delta_BS_RIS(2)^2));
    angles.BS_AoD.distance = norm(delta_BS_RIS);

    %% ============================
    % 2. RIS AoA from BS
    % ============================
    angles.RIS_AoA.azimuth = mod(angles.BS_AoD.azimuth + 180, 360);
    angles.RIS_AoA.elevation = -angles.BS_AoD.elevation;
    angles.RIS_AoA.distance = angles.BS_AoD.distance;

    %% ============================
    % 3. RIS -> users/BD angles
    % ============================
    for c = 1:Nclusters
        for k = 1:Kusers

            user = users{c}(k,:);
            delta_RIS_user = user - RIS;

            azimuth = atan2d(delta_RIS_user(2), delta_RIS_user(1));
            elevation = atan2d(delta_RIS_user(3), ...
                sqrt(delta_RIS_user(1)^2 + delta_RIS_user(2)^2));

            angles.RIS_AoD.clusters(c).users(k).user = sprintf('C%d-U%d', c, k);
            angles.RIS_AoD.clusters(c).users(k).azimuth = azimuth;
            angles.RIS_AoD.clusters(c).users(k).elevation = elevation;
            angles.RIS_AoD.clusters(c).users(k).distance = norm(delta_RIS_user);
            angles.RIS_AoD.clusters(c).users(k).position = user;
        end
    end

    %% ============================
    % 4. Distances
    % ============================
    distances = struct();

    for c = 1:Nclusters

        BD = users{c}(3,:);
        clusterDistances = struct();

        for u = 1:2
            user_pos = users{c}(u,:);
            clusterDistances.(['User' num2str(u) '_BD']) = norm(user_pos - BD);
            clusterDistances.(['RIS_to_User' num2str(u)]) = norm(user_pos - RIS);
            clusterDistances.(['BS_to_User' num2str(u)]) = norm(user_pos - BS);
        end

        clusterDistances.User1_User2 = norm(users{c}(1,:) - users{c}(2,:));
        clusterDistances.RIS_to_BD = norm(BD - RIS);
        clusterDistances.BS_to_BD = norm(BD - BS);

        distances.(['Cluster' num2str(c)]) = clusterDistances;
    end

    %% ============================
    % 5. Visualization
    % ============================
    figure; hold on; grid on; axis equal;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('Short-Distance BS, RIS, Two NOMA Clusters and BDs');

    scatter3(BS(1),BS(2),BS(3),100,'r','filled','^');
    text(BS(1),BS(2),BS(3)+0.5,'BS');

    scatter3(RIS(1),RIS(2),RIS(3),100,'b','filled','s');
    text(RIS(1),RIS(2),RIS(3)+0.5,'RIS');

    colors = lines(Nclusters);

    for c = 1:Nclusters

        center = cluster_centers(c,:);
        scatter3(center(1),center(2),center(3),80,colors(c,:),'d','filled');
        text(center(1),center(2),center(3)+0.5,sprintf('Cluster %d',c));

        theta = linspace(0,2*pi,100);
        x_circle = center(1) + radii(c)*cos(theta);
        y_circle = center(2) + radii(c)*sin(theta);
        z_circle = center(3)*ones(size(theta));
        plot3(x_circle,y_circle,z_circle,'--','Color',colors(c,:),'LineWidth',1.2);

        cluster_users = users{c};

        % NOMA users
        scatter3(cluster_users(1:2,1),cluster_users(1:2,2),cluster_users(1:2,3), ...
            60,colors(c,:),'filled','o');

        % BD
        scatter3(cluster_users(3,1),cluster_users(3,2),cluster_users(3,3), ...
            90,colors(c,:),'filled','diamond');

        text(cluster_users(1,1),cluster_users(1,2),cluster_users(1,3)+0.5,'U1');
        text(cluster_users(2,1),cluster_users(2,2),cluster_users(2,3)+0.5,'U2');
        text(cluster_users(3,1),cluster_users(3,2),cluster_users(3,3)+0.5,'BD');

    end

    view(45,25);
    hold off;

    %% Print distances
    fprintf('\n=== Distance Summary ===\n');
    fprintf('BS-RIS distance: %.2f m\n', norm(RIS - BS));

    for c = 1:Nclusters
        fprintf('\nCluster %d:\n', c);
        fprintf('RIS-U1: %.2f m\n', distances.(['Cluster' num2str(c)]).RIS_to_User1);
        fprintf('RIS-U2: %.2f m\n', distances.(['Cluster' num2str(c)]).RIS_to_User2);
        fprintf('RIS-BD: %.2f m\n', distances.(['Cluster' num2str(c)]).RIS_to_BD);
        fprintf('U1-BD: %.2f m\n', distances.(['Cluster' num2str(c)]).User1_BD);
        fprintf('U2-BD: %.2f m\n', distances.(['Cluster' num2str(c)]).User2_BD);
    end

end