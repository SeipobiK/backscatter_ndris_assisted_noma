function [values] = para_init()
    % Author: Kgomotjo Seipobi

    % ===============================
    % Basic Parameter2
    % ===============================
    values.noise_dB = -90; % noise power spectral density in dBW/Hz
    values.scall=400;
    % values.noise = 10^((val2es.noise_dB - 30)/10) * (62500000000); 
    % values.noise = 10^(values.noise_dB/10)* (200)^4; 
    values.noise = (1e-12) * (400)^4; % override with fixed scaling

    values.alpha_k_n = 0.5;
    values.alpha_k_f = 0.5;

    % User weights
    values.weights_n = 1; % near user
    values.weights_f = 1; % far user
    values.weights_c = 1;    % backscatter
    values.weights = [0.2, 0.3, 0.3, 0.2]; % example custom weights

    % System dimensions
    values.K_u = 5; % users per cluster
    values.K_c = 2; % NOMA users per cluster
    values.K = 1;   % number of clusters
    values.M =8;   % antennas at BS
    values.RIS_size = [2,16]; % RIS dimensions (rows, columns)
    values.N = prod(values.RIS_size);  

    % Power and rate requirements

    values.P_max =100; % maximum transmit power in Watts
    values.eta = 0.7; % backscatter coefficient
    values.R_min_f = 0.01; 
    values.R_min_n = 0.01; 
    values.R_c_min = 0.01; 
    values.R_min_c = 0.01;
    % values.R_c_min = 0.1; 
    values.nu_n = 1; 


    values.nu_f = 1; 
    values.nu_c = 1; 
    

    % Iterations
    values.MC_MAX =10;
    values.outer_iter = 15; 
    values.max_iter = 10; 
    values.tol = 1e-5; 

    values.pathloss = @(d) 30 + 22*log10(d); 
    values.rician = 10^(-10/10); 

    % ===============================
    % Geometry: BS, RIS, Users
    % ===============================
    [angles, users,distances] = calculate_angles();   

    values.BSTdist=distances;
    % BS location (AoD)
    values.BS_loc = [angles.BS_AoD.distance, ...
                     angles.BS_AoD.elevation, ...
                     angles.BS_AoD.azimuth];

    % RIS location (AoA)
    values.RIS_loc = [angles.RIS_AoA.distance, ...
                      angles.RIS_AoA.elevation, ...
                      angles.RIS_AoA.azimuth];

    % Store users per cluster
    values.users = users;

    % ===============================
    % User locations (d, elev, az) wrt RIS
    % ===============================
    values.userloc = zeros(values.K_u, values.K, 3);  
    for c = 1:values.K
        for u = 1:values.K_u
            values.userloc(u, c, 1) = angles.RIS_AoD.clusters(c).users(u).distance; 
            values.userloc(u, c, 2) = angles.RIS_AoD.clusters(c).users(u).elevation; 
            values.userloc(u, c, 3) = angles.RIS_AoD.clusters(c).users(u).azimuth;   
        end
    end
end

