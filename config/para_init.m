function [values] = para_init()
    % Author: Kgomotjo Seipobi

    % ===============================
    % Basic Parameters
    % ===============================
    values.noise_dB = -90;
    values.scall    = 500;
    values.noise    = (1e-12) * (500)^4;

    values.alpha_k_n = 0.55;
    values.alpha_k_f = 0.45;
    values.FT        = 2;
    values.rho       = 0.6;
    values.epsilon   = 1e-6;
    values.bst_threshold = 0.000010000;

    % User weights
    values.weights_n = 1;
    values.weights_f = 1;
    values.weights_c = 1;
    values.weights   = [0.2, 0.3, 0.3, 0.2];

    % System dimensions
    values.K_u      = 3;
    values.K_c      = 2;
    values.K        = 2;
    values.M        = 8;
    values.RIS_size = [2, 16];
    values.N        = prod(values.RIS_size);

    % Power and rate requirements
    values.P_max  = 10;
    values.eta    = [0.5, 0.5];
    values.R_min_n = 0.5;
    values.R_c_min = 0.5;
    values.R_min_c = 0.5;
    values.nu_n   = 3;
    values.nu_f   = 1;
    values.nu_c   = 1;

    % WSR weights (equal weighting)
    values.omega_user = ones(values.K, values.K_c);
    values.omega_back = ones(values.K, 1);

    % Iterations (your existing fields)
    values.MC_MAX    = 1000;
    values.outer_iter = 20;
    values.max_iter  = 15;
    values.tol       = 1e-5;

    % ===============================
    % Passive Beamforming Parameters
    % (passive_beamforming_riemannian)
    % ===============================

    % — mapped from your existing fields —
    % outer_iter already exists  →  used as max_outer
    % max_iter   already exists  →  used as max_inner
    % tol        already exists  →  used as tol_inner
    % phi_init   already exists  →  used as phi_init
    % These aliases keep the rest of your code unchanged:
    values.max_outer  = values.outer_iter;   % 5  outer dual updates
    values.max_inner  = values.max_iter;     % 15 Riemannian CG steps
    values.tol_outer  = 1e-4;               % relative L change (outer)
    values.tol_inner  = values.tol;         % ||G_M||_F threshold (inner)
    values.phi_init   =1e-2;               % unit step (standard for manifold CG)

    % — new fields —
    values.phi_min    = 1e-8;   % give up line search below this step size
    values.c_armijo   = 1e-4;   % sufficient increase constant (Nocedal & Wright)

    % Symmetry penalty for BD-RIS reciprocity  Theta = Theta^T
    values.nu_init    = 0;   % start much smaller
    values.nu_max     = 0;   % much lower cap
    values.rho_nu     = 1.2;    % grow slowly, not doubling
    values.delta_nu   = 0.5;    % easier reduction threshold

    % Dual variable step schedule: beta_t = beta_0_dual / sqrt(t)
    values.beta_0_dual = 0;   % gives ~0.014 at iteration 50

    % ===============================
    % Path loss and channel
    % ===============================
    values.pathloss = @(d) 30 + 22*log10(d);
    values.rician   = 10^(-10/10);

    % ===============================
    % Geometry: BS, RIS, Users
    % ===============================
    [angles, users, distances] = calculate_angles();

    values.BSTdist  = distances;

    values.BS_loc   = [angles.BS_AoD.distance, ...
                       angles.BS_AoD.elevation, ...
                       angles.BS_AoD.azimuth];

    values.RIS_loc  = [angles.RIS_AoA.distance, ...
                       angles.RIS_AoA.elevation, ...
                       angles.RIS_AoA.azimuth];

    values.users    = users;

    values.userloc  = zeros(values.K_u, values.K, 3);
    for c = 1:values.K
        for u = 1:values.K_u
            values.userloc(u,c,1) = angles.RIS_AoD.clusters(c).users(u).distance;
            values.userloc(u,c,2) = angles.RIS_AoD.clusters(c).users(u).elevation;
            values.userloc(u,c,3) = angles.RIS_AoD.clusters(c).users(u).azimuth;
        end
    end
end