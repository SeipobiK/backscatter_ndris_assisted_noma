function [Theta_opt, nu_opt, duals_vars, duals_vars_c] = ...
    passive_beamforming_riemannian(para, w_k, Theta_init, channel_data, ...
                                   decoding_order, alpha, eta)
% PASSIVE_BEAMFORMING_RIEMANNIAN
%   Riemannian conjugate gradient ascent for the BD-RIS scattering matrix.
%
%   Implements Algorithm (P_pass):
%     - Outer loop : dual variable updates (mu_{k,i}, mu_{k,b}, nu)
%     - Inner loop : Riemannian CGA on the unitary manifold U(N)
%                   with QR retraction and Polak-Ribiere parameter
%
%   Notation matches dissertation eqs:
%     (passive_total_gradient)   — Euclidean gradient
%     (riemannian_projection)    — tangent-space projection
%     (ascent_direction)         — CG ascent direction
%     (polak_ribiere)            — PR beta with max(0,.) clamp
%     (retraction)               — QR retraction
%     (vector_transport)         — projection-based transport
%     (nu_update)                — symmetry penalty weight update
%     (passive_mu_user_update)   — NOMA QoS dual update
%     (passive_mu_backscatter_update) — backscatter QoS dual update
%
% INPUTS
%   para          : struct with fields
%                     .K, .K_c          — cluster / user counts
%                     .R_min_n          — min NOMA rate (bits/s/Hz)
%                     .R_min_b          — min backscatter rate
%                     .noise            — noise variance sigma^2
%                     .nu_init          — initial symmetry penalty weight
%                     .nu_max           — cap on nu
%                     .rho_nu           — nu growth factor  (e.g. 2)
%                     .delta_nu         — required symmetry reduction ratio (e.g. 0.25)
%                     .omega_user       — K x K_c WSR weights (set to ones if equal)
%                     .omega_back       — K x 1   WSR weights
%                     .max_outer        — outer iteration limit
%                     .max_inner        — inner CG iteration limit
%                     .tol_outer        — outer convergence tolerance
%                     .tol_inner        — inner CG convergence tolerance
%                     .c_armijo         — Armijo constant (e.g. 1e-4)
%                     .phi_init         — initial step size (e.g. 1)
%                     .phi_min          — minimum step size before giving up
%   w_k           : N_t x K  beamforming matrix (fixed during this subproblem)
%   Theta_init    : N  x N   initial scattering matrix (unitary)
%   channel_data  : struct with fields
%                     .H_all   N x N_t  composite BS-RIS channel
%                     .g{k,i}  N x 1    RIS-to-user
%                     .g_b{k}  N x 1    RIS-to-backscatter-decoder
%                     .f{k,i}  scalar   tag reflection coefficient
%   decoding_order : K x K_c  SIC order
%   alpha          : K x K_c  power allocation coefficients
%   eta            : K x 1    backscatter reflection efficiencies
%
% OUTPUTS
%   Theta_opt     : N x N   optimised scattering matrix
%   nu_opt        : scalar  final symmetry penalty weight
%   duals_vars    : K x K_c final NOMA QoS multipliers
%   duals_vars_c  : K x 1  final backscatter QoS multipliers

%% -----------------------------------------------------------------------
%  Unpack and initialise
%% -----------------------------------------------------------------------
K   = para.K;
K_c = para.K_c;
N   = size(Theta_init, 1);

% Symmetry penalty
nu          =  para.nu_init;
sym_viol_prev = Inf;           % ||Theta - Theta^T||_F from previous outer iter

% Dual variables — initialise to zero
duals_vars   = zeros(K, K_c);
duals_vars_c = zeros(K, 1);

% Diminishing step schedule for dual updates: beta_t = beta_0 / sqrt(t)
beta_0_dual  = 1e-3;

Theta = Theta_init;

%% -----------------------------------------------------------------------
%  OUTER LOOP — dual variable updates



%% -----------------------------------------------------------------------
for t_out = 1 : para.max_iter

      J_r=eye(N); J_t=eye(N);

      [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);

                             [sum_rate, R, R_c, noma_signal, BST_signal, noma_interference, BST_interference, ...
            intra_i, inteer_i, inteer_b, inteer_b_all] = compute_wsr(...
            para, H, H_c, decoding_order, alpha, w_k,eta);

            % At the very start of the outer loop, before anything else:
        fprintf('=== Outer iter %d ===\n', t_out);
        fprintf('  WSR before inner CG: %.4f\n', ...
            sum_rate);





    %% -------------------------------------------------------------------
    %  Step 1 — Recompute cascaded channels with current Theta
    %            H{k,i}   = g_{k,i}^H Theta H   (1 x N_t row vector)
    %            H_c{k,i} = f_{k,i} g_{k,b}^H Theta H  (1 x N_t row vector)
    %  These are fed into compute_fp_dual_vars which expects them pre-formed.
    %% -------------------------------------------------------------------
    H_cas  = cell(K, K_c);   % cascaded NOMA channels
    H_c    = cell(K, K_c);   % cascaded backscatter channels

    for k = 1:K
        for i = 1:K_c
            H_cas{k,i} = channel_data.g{k,i}' * Theta * channel_data.H_all;
            H_c{k,i}   = channel_data.f{k,i}  * ...
                          channel_data.g_b{k}' * Theta * channel_data.H_all;
        end
    end

    %% -------------------------------------------------------------------
    %  Step 2 — Update FP auxiliary vars (lambda, y, z) and dual vars
    %            using current Theta's cascaded channels
    %% -------------------------------------------------------------------
    learning_rate_t = beta_0_dual / sqrt(t_out);   % diminishing dual step

    [beta_fp, zeta_fp, y_fp, z_fp, ~, ~, duals_vars, duals_vars_c] = ...
        compute_fp_gual_vars(para, H_cas, H_c, decoding_order, alpha, w_k, eta, ...
                             duals_vars, learning_rate_t, duals_vars_c);

    %% -------------------------------------------------------------------
    %  Step 3 — Inner Riemannian CGA to update Theta
    %% -------------------------------------------------------------------
    Delta_prev  = zeros(N, N);     % previous CG direction (tangent vector)
    G_prev      = zeros(N, N);     % previous Riemannian gradient
    obj_prev    = -Inf;

    for t_in = 1 : para.max_iter

              [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);

                             [sum_rate, R, R_c, noma_signal, BST_signal, noma_interference, BST_interference, ...
            intra_i, inteer_i, inteer_b, inteer_b_all] = compute_wsr(...
            para, H, H_c, decoding_order, alpha, w_k,eta);

            aug_val   = augmented_lagrangian(para, w_k, Theta, channel_data, ...
                decoding_order, beta_fp, zeta_fp, y_fp, z_fp, ...
                duals_vars, duals_vars_c, alpha, eta, nu);
            fprintf('    inner %d: WSR=%.4f  AugLag=%.6f\n', t_in, sum_rate, aug_val);
            %% Refresh cascaded channels at current Theta
            for k = 1:K
                for i = 1:K_c
                    H_cas{k,i} = channel_data.g{k,i}' * Theta * channel_data.H_all;
                    H_c{k,i}   = channel_data.f{k,i} * ...
                                channel_data.g_b{k}' * Theta * channel_data.H_all;
                end
            end

                %% Refresh FP auxiliaries only — dual vars held fixed (lr=0)
                [beta_fp, zeta_fp, y_fp, z_fp, ~, ~, ~, ~] = ...
                    compute_fp_gual_vars(para, H_cas, H_c, decoding_order, alpha, w_k, eta, ...
                                        duals_vars, 0, duals_vars_c);
            % TEMP DIAGNOSTIC — paste right after auxiliaries are refreshed, before gradient
            for k = 1:K
                for i = 1:K_c
                    lam  = beta_fp(k,i);
                    y_ki = y_fp(k,i);
                    
                    % recompute A_ki, B_ki exactly as augmented_lagrangian does
                    g_ki = channel_data.g{k,i};
                    d_ki = g_ki' * Theta * channel_data.H_all * w_k(:,k);
                    inter = 0;
                    for j = 1:K
                        if j ~= k
                            inter = inter + abs(g_ki' * Theta * channel_data.H_all * w_k(:,j))^2;
                        end
                    end
                    order_k = decoding_order(k,:);
                    pos_i = find(order_k == i);
                    alpha_sum = 0;
                    for idx = pos_i+1:K_c
                        alpha_sum = alpha_sum + alpha(k, order_k(idx));
                    end
                    intra = abs(d_ki)^2 * alpha_sum;
                    bs_int = 0;
                    for kp = 1:K
                        f_kpi = channel_data.f{kp,i};
                        g_kpb = channel_data.g_b{kp};
                        d_kpbi = f_kpi * (g_kpb' * Theta * channel_data.H_all * w_k(:,kp));
                        bs_int = bs_int + eta(kp) * abs(d_kpbi)^2;
                    end
                    A_ki = abs(d_ki)^2 * alpha(k,i);
                    B_ki = inter + intra + bs_int + para.noise;
                    
                    term1 = log(1+lam) - lam;
                    term2 = 2*real(conj(y_ki)*sqrt((1+lam)*alpha(k,i))*d_ki);
                    term3 = -abs(y_ki)^2 * (A_ki + B_ki);
                    
                    fprintf('k=%d i=%d: lam=%.4e y=%.4e A=%.4e B=%.4e | t1=%.4e t2=%.4e t3=%.4e\n', ...
                        k, i, lam, abs(y_ki), A_ki, B_ki, term1, term2, term3);
                end
            end

        %% ---------------------------------------------------------------
        %  3a. Euclidean gradient  (eq. passive_total_gradient)
        %% ---------------------------------------------------------------
        E_grad = compute_eucl_grad(para, w_k, Theta, channel_data, ...
                                   decoding_order, beta_fp, zeta_fp, ...
                                   y_fp, z_fp, duals_vars, duals_vars_c, ...
                                   alpha, eta);

        %% ---------------------------------------------------------------
        %  3b. Riemannian gradient via tangent-space projection
        %      (eq. riemannian_projection)
        %
        %      G_M = E_grad - Theta * sym( Theta^H * E_grad )
        %      sym(X) = (X + X^H) / 2
        %% ---------------------------------------------------------------
        G_M = riemannian_grad(Theta, E_grad);

        %% ---------------------------------------------------------------
        %  3c. Polak-Ribiere parameter  (eq. polak_ribiere)
        %
        %      epsilon = max(0,  <G_M, G_M - T(G_prev)> / ||G_prev||_F^2 )
        %
        %  where T(G_prev) = transport of previous gradient to current
        %  tangent space = Pi_{Theta}(G_prev)
        %% ---------------------------------------------------------------
        if t_in == 1
            epsilon_pr = 0;
        else
            T_G_prev      = proj_tangent(Theta, G_prev);   % transport G_prev
            numer         = real(trace(G_M' * (G_M - T_G_prev)));
            denom         = real(trace(G_prev' * G_prev));
            epsilon_pr    = max(0, numer / (denom + eps));
        end

        %% ---------------------------------------------------------------
        %  3d. CG ascent direction  (eq. ascent_direction)
        %
        %      Delta = +G_M + epsilon * T(Delta_prev)
        %% ---------------------------------------------------------------
        if t_in == 1
            Delta = G_M;
        else
            T_Delta_prev = proj_tangent(Theta, Delta_prev);
            Delta        = G_M + epsilon_pr * T_Delta_prev;
        end

        %  3e. Armijo backtracking line search
        %% ---------------------------------------------------------------
        obj_cur = augmented_lagrangian(para, w_k, Theta, channel_data, ...
                                    decoding_order, beta_fp, zeta_fp, ...
                                    y_fp, z_fp, duals_vars, duals_vars_c, ...
                                    alpha, eta, nu);

        slope = real(trace(G_M' * Delta));   % Riemannian inner product
        if slope <= 0
            Delta = G_M;                     % reset to steepest ascent
            slope = real(trace(G_M' * G_M));
        end

        phi       = para.phi_init;
        armijo_ok = false;

        while phi >= para.phi_min
            Theta_new = polar_retract(Theta, phi * Delta);
            
            % Recompute auxiliaries at Theta_new for honest evaluation
            for k = 1:K
                for i = 1:K_c
                    H_cas_new{k,i} = channel_data.g{k,i}' * Theta_new * channel_data.H_all;
                    H_c_new{k,i}   = channel_data.f{k,i} * ...
                                    channel_data.g_b{k}' * Theta_new * channel_data.H_all;
                end
            end
            [beta_new, zeta_new, y_new, z_new, ~, ~, ~, ~] = ...
                compute_fp_gual_vars(para, H_cas_new, H_c_new, decoding_order, ...
                                    alpha, w_k, eta, duals_vars, 0, duals_vars_c);
            
            obj_new = augmented_lagrangian(para, w_k, Theta_new, channel_data, ...
                        decoding_order, beta_new, zeta_new, y_new, z_new, ...
                        duals_vars, duals_vars_c, alpha, eta, 0);
            
            if obj_new >= obj_cur + para.c_armijo * phi * slope
                beta_fp  = beta_new;   % accept new auxiliaries too
                zeta_fp  = zeta_new;
                y_fp     = y_new;
                z_fp     = z_new;
                armijo_ok = true;
                break;
            end
            phi = phi / 2;
        end

        if ~armijo_ok
            % Line search failed — fall back to gradient step with min step
            Theta_new = polar_retract(Theta, para.phi_min * G_M);
        end

        %% ---------------------------------------------------------------
        %  3f. Accept step
        %% ---------------------------------------------------------------
        G_prev     = G_M;
        Delta_prev = Delta;
        Theta      = Theta_new;

        %% ---------------------------------------------------------------
        %  3g. Inner convergence check
        %% ---------------------------------------------------------------
        grad_norm = norm(G_M, 'fro');
        if grad_norm < para.tol_inner
            break;
        end

    end  % inner CG loop

    %% -------------------------------------------------------------------
    %  Step 4 — Update symmetry penalty weight nu  (eq. nu_update)
    %% -------------------------------------------------------------------
    sym_viol_cur = norm(Theta - Theta.', 'fro');

    if sym_viol_cur > para.delta_nu * sym_viol_prev
        nu = min(para.rho_nu * nu, para.nu_max);
    end
    sym_viol_prev = sym_viol_cur;

    %% -------------------------------------------------------------------
    %  Step 5 — Outer convergence check
    %% -------------------------------------------------------------------
    obj_check = augmented_lagrangian(para, w_k, Theta, channel_data, ...
                                     decoding_order, beta_fp, zeta_fp, ...
                                     y_fp, z_fp, duals_vars, duals_vars_c, ...
                                     alpha, eta, nu);
    if abs(obj_check - obj_prev) / (abs(obj_prev) + eps) < para.tol_outer
        break;
    end


    obj_prev = obj_check;

          [H, H_c] = build_effective_channels(para, channel_data, Theta, J_r, J_t);

                             [sum_rate, R, R_c, noma_signal, BST_signal, noma_interference, BST_interference, ...
            intra_i, inteer_i, inteer_b, inteer_b_all] = compute_wsr(...
            para, H, H_c, decoding_order, alpha, w_k,eta);

            % At the very start of the outer loop, before anything else:
        fprintf('=== Outer iter %d ===\n', t_out);
        fprintf('  WSR before inner CG: %.4f\n', ...
            sum_rate);

            % At the end of the outer loop, after inner CG:
        fprintf('  WSR after inner CG:  %.4f\n', ...
            sum_rate);
        fprintf('  norm(duals_vars):    %.4e\n', norm(duals_vars,'fro'));
        fprintf('  norm(duals_vars_c):  %.4e\n', norm(duals_vars_c));
        fprintf('  sym_viol:            %.4e\n', norm(Theta - Theta.','fro'));
        fprintf('  nu:                  %.4e\n', nu);

end  % outer loop

Theta_opt = Theta;
nu_opt    = nu;

end  % main function


%% =======================================================================
%  LOCAL HELPER FUNCTIONS
%% =======================================================================

function G_M = riemannian_grad(Theta, E_grad)
% RIEMANNIAN_GRAD  Project Euclidean gradient onto tangent space of U(N).
%   (eq. riemannian_projection)
%   G_M = E_grad - Theta * sym( Theta^H * E_grad )
%   sym(X) = (X + X^H) / 2
    S   = Theta' * E_grad;
    G_M = E_grad - Theta * ((S + S') / 2);
end


function P = proj_tangent(Theta, Xi)
% PROJ_TANGENT  Project ambient matrix Xi onto tangent space at Theta.
%   Used for both vector transport and Riemannian gradient.
%   (eq. vector_transport)
%   P = Xi - Theta * sym( Theta^H * Xi )
    S = Theta' * Xi;
    P = Xi - Theta * ((S + S') / 2);
end


function Theta_new = polar_retract(Theta, V)
% POLAR_RETRACT  Retract tangent update onto unitary manifold via
%                polar decomposition.  More stable than QR for square
%                matrices — no sign ambiguity, nearest unitary in
%                Frobenius norm.
%   Theta_new = U * Vh  where  Theta + V = U * S * Vh  (thin SVD)
    [U, ~, Vh] = svd(Theta + V, 'econ');
    Theta_new  = U * Vh;
end


function val = augmented_lagrangian(para, w_k, Theta, channel_data, ...
                                     decoding_order, beta_fp, zeta_fp, ...
                                     y_fp, z_fp, duals_vars, duals_vars_c, ...
                                     alpha, eta, nu)
% AUGMENTED_LAGRANGIAN  Evaluate L_pass at current Theta.
%   (eq. passive_aug_lag)
%   L = f_QT + sum mu_{k,i}*(A_{k,i} - gbar*B_{k,i})
%             + sum mu_{k,b}*(A_{k,b} - gbar*B_{k,b})
%             - nu*||Theta - Theta^T||_F^2

    K   = para.K;
    K_c = para.K_c;
    H   = channel_data.H_all;
    ln2 = log(2);

    gamma_noma = 2^para.R_min_n - 1;
    gamma_back = 2^para.R_c_min - 1;

    val = 0;

    for k = 1:K
        order_k     = decoding_order(k,:);
        strong_user = order_k(end);

        for i = 1:K_c
            g_ki = channel_data.g{k,i};

            % Cascaded channels
            d_ki = g_ki' * Theta * H * w_k(:,k);   % desired

            % Inter-cluster interference
            inter = 0;
            for j = 1:K
                if j ~= k
                    inter = inter + abs(g_ki' * Theta * H * w_k(:,j))^2;
                end
            end

            % Intra-cluster interference
            pos_i     = find(order_k == i);
            alpha_sum = 0;
            for idx = pos_i+1 : K_c
                alpha_sum = alpha_sum + alpha(k, order_k(idx));
            end
            intra = abs(d_ki)^2 * alpha_sum;

            % Backscatter interference on user (k,i)
            bs_int = 0;
            for kp = 1:K
                f_kpi  = channel_data.f{kp,i};
                g_kpb  = channel_data.g_b{kp};
                d_kpbi = f_kpi * (g_kpb' * Theta * H * w_k(:,kp));
                bs_int = bs_int + eta(kp) * abs(d_kpbi)^2;
            end

            A_ki = abs(d_ki)^2 * alpha(k,i);
            B_ki = inter + intra + bs_int + para.noise;

            % QT user term  (omega / ln2) * [ ln(1+lambda) - lambda
            %                + 2Re{y* sqrt((1+lambda)*alpha) d} - |y|^2*(A+B) ]
            lam  = beta_fp(k,i);
            y_ki = y_fp(k,i);
            w_u  = para.omega_user(k,i) / ln2;
            val  = val + w_u * ( log(1+lam) - lam ...
                       + 2*real(conj(y_ki)*sqrt((1+lam)*alpha(k,i))*d_ki) ...
                       - abs(y_ki)^2 * (A_ki + B_ki) );

            % NOMA QoS dual term
            val = val + duals_vars(k,i) * (A_ki - gamma_noma * B_ki);

            %% Backscatter decoder (strong user only)
            if i == strong_user
                f_ki_star = channel_data.f{k,i};
                g_kb      = channel_data.g_b{k};
                d_kb      = f_ki_star * (g_kb' * Theta * H * w_k(:,k));

                % B2B + inter interference for backscatter decoder
                inter_bc = 0;
                for j = 1:K
                    if j ~= k
                        d_kb_j   = f_ki_star * (g_kb' * Theta * H * w_k(:,j));
                        inter_bc = inter_bc + abs(d_kb_j)^2;
                    end
                end
                b2b = 0;
                for kp = 1:K
                    if kp ~= k
                        f_kp  = channel_data.f{kp, decoding_order(kp,end)};
                        g_kpb = channel_data.g_b{kp};
                        d_kpb = f_kp * (g_kpb' * Theta * H * w_k(:,kp));
                        b2b   = b2b + eta(kp) * abs(d_kpb)^2;
                    end
                end

                A_kb = eta(k) * abs(d_kb)^2;
                B_kb = inter_bc + b2b + para.noise;

                lam_b = zeta_fp(k);
                z_k   = z_fp(k);
                w_b   = para.omega_back(k) / ln2;
                val   = val + w_b * ( log(1+lam_b) - lam_b ...
                            + 2*real(conj(z_k)*sqrt((1+lam_b)*eta(k))*d_kb) ...
                            - abs(z_k)^2 * (A_kb + B_kb) );

                % Backscatter QoS dual term
                val = val + duals_vars_c(k) * (A_kb - gamma_back * B_kb);
            end
        end
    end

    % Symmetry penalty  -nu * ||Theta - Theta^T||_F^2
    % val = val - nu * norm(Theta - Theta.', 'fro')^2;

end