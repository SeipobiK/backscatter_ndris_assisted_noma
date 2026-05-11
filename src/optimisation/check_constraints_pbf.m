function [V_opt, A_opt, B_opt, A_c_opt, B_c_opt, obj_prev, status, slack_report, slacks] = ...
    check_constraints_pbf(para, w_k, channel_data, decoding_order, ...
    A_prev, B_prev, A_c_prev, B_c_prev, alpha, J_r, J_t, eta)

    %% Parameters
    K       = para.K;
    K_c     = para.K_c;
    N       = para.N;
    noise   = para.noise;
    R_min   = para.R_min_n;
    R_min_c = para.R_min_c;
    eh      = para.bst_threshold;
    rho     = para.rho;

    %% Precompute effective channel covariance matrices
    H   = cell(K, K_c);
    H_c = cell(K, K_c);

    for kk = 1:K
        for uu = 1:K_c

            h_eff = diag(channel_data.g{kk,uu}' * J_r) * ...
                    J_t * channel_data.H_all * w_k(:,kk);

            h_c_eff = diag(channel_data.g_b{kk}' * J_r) * ...
                      J_t * channel_data.H_all * ...
                      channel_data.f{kk,uu} * w_k(:,kk);

            H{kk,uu}   = h_eff * h_eff';
            H_c{kk,uu} = h_c_eff * h_c_eff';
        end
    end

    %% CVX problem
    cvx_begin quiet
        cvx_solver mosek

        %% Optimization variables
        variable V(N,N) hermitian semidefinite

        variable A(K,K_c) nonnegative
        variable B(K,K_c) nonnegative
        variable A_c(K) nonnegative
        variable B_c(K) nonnegative

        variable R(K,K_c) nonnegative
        variable R_c(K) nonnegative

        %% Slack variables
        variable s_diag(N) nonnegative

        variable s_rate(K,K_c) nonnegative
        variable s_qos(K,K_c) nonnegative
        variable s_sig(K,K_c) nonnegative
        variable s_int(K,K_c) nonnegative

        variable s_eh(K) nonnegative
        variable s_rate_c(K) nonnegative
        variable s_sig_c(K) nonnegative
        variable s_int_c(K) nonnegative
        variable s_qos_c(K) nonnegative

        %% Objective: minimise total constraint violation
        minimize( ...
            sum(s_diag(:)) + ...
            sum(s_rate(:)) + ...
            sum(s_qos(:)) + ...
            sum(s_sig(:)) + ...
            sum(s_int(:)) + ...
            sum(s_eh(:)) + ...
            sum(s_rate_c(:)) + ...
            sum(s_sig_c(:)) + ...
            sum(s_int_c(:)) + ...
            sum(s_qos_c(:)) )

        subject to

            %% Relaxed unit-modulus diagonal constraints
            for mm = 1:N
                V(mm,mm) <= 1 + s_diag(mm);
                V(mm,mm) >= 1 - s_diag(mm);
            end

            %% Main constraints
            for kk = 1:K

                order_k = decoding_order(kk,:);
                strong_user = order_k(end);

                for uu = 1:K_c

                    %% NOMA SCA rate approximation
                    rate_sca = log2(1 + 1/(A_prev(kk,uu) * B_prev(kk,uu))) ...
                        - (log2(exp(1)) / ...
                        (A_prev(kk,uu) * (1 + A_prev(kk,uu) * B_prev(kk,uu)))) ...
                        * (A(kk,uu) - A_prev(kk,uu)) ...
                        - (log2(exp(1)) / ...
                        (B_prev(kk,uu) * (1 + A_prev(kk,uu) * B_prev(kk,uu)))) ...
                        * (B(kk,uu) - B_prev(kk,uu));

                    R(kk,uu) <= rate_sca + s_rate(kk,uu);

                    %% NOMA QoS
                    R(kk,uu) + s_qos(kk,uu) >= R_min;

                    %% Inter-cluster interference
                    inter = 0;
                    inter_b = 0;

                    for jj = 1:K
                        if jj ~= kk

                            h_inter = diag(channel_data.g{kk,uu}' * J_r) * ...
                                      J_t * channel_data.H_all * w_k(:,jj);

                            inter = inter + real(trace(V * (h_inter * h_inter')));

                            h_inter_b = diag(channel_data.g_b{kk}' * J_r) * ...
                                        J_t * channel_data.H_all * ...
                                        channel_data.f{kk,uu} * w_k(:,jj);

                            inter_b = inter_b + ...
                                eta(jj) * real(trace(V * (h_inter_b * h_inter_b')));
                        end
                    end

                    %% Intra-cluster NOMA interference
                    intra = 0;
                    pos_user = find(order_k == uu);

                    for order_pos = pos_user+1:K_c
                        interfering_user = order_k(order_pos);
                        intra = intra + alpha(kk,interfering_user) * ...
                                real(trace(V * H{kk,uu}));
                    end

                    %% NOMA signal constraint
                    signal_noma = alpha(kk,uu) * real(trace(V * H{kk,uu}));

                    inv_pos(A(kk,uu)) <= signal_noma + s_sig(kk,uu);

                    %% NOMA interference constraint
                    backscatter_self_int = eta(kk) * real(trace(V * H_c{kk,uu}));

                    B(kk,uu) + s_int(kk,uu) >= ...
                        intra + inter + inter_b + backscatter_self_int + noise;

                    %% Backscatter constraints only at strong/SIC user
                    if uu == strong_user

                        %% Energy harvesting constraint
                        h_eh = diag(channel_data.g_b{kk}' * J_r) * ...
                               J_t * channel_data.H_all * w_k(:,kk);

                        harvested_energy = (1 - eta(kk)) * rho * ...
                            real(trace(V * (h_eh * h_eh')));

                        harvested_energy + s_eh(kk) >= eh;

                        %% Backscatter signal constraint
                        signal_c = eta(kk) * real(trace(V * H_c{kk,uu}));

                        inv_pos(A_c(kk)) <= signal_c + s_sig_c(kk);

                        %% Backscatter SCA rate approximation
                        rate_c_sca = log2(1 + 1/(A_c_prev(kk) * B_c_prev(kk))) ...
                            - (log2(exp(1)) / ...
                            (A_c_prev(kk) * (1 + A_c_prev(kk) * B_c_prev(kk)))) ...
                            * (A_c(kk) - A_c_prev(kk)) ...
                            - (log2(exp(1)) / ...
                            (B_c_prev(kk) * (1 + A_c_prev(kk) * B_c_prev(kk)))) ...
                            * (B_c(kk) - B_c_prev(kk));

                        R_c(kk) <= rate_c_sca + s_rate_c(kk);

                        %% Backscatter interference constraint
                        B_c(kk) + s_int_c(kk) >= inter + inter_b + noise;

                        %% Backscatter QoS
                        R_c(kk) + s_qos_c(kk) >= R_min_c;
                    end
                end
            end

    cvx_end

    %% Outputs
    obj_prev = cvx_optval;
    status   = cvx_status;

    V_opt   = V;
    A_opt   = A;
    B_opt   = B;
    A_c_opt = A_c;
    B_c_opt = B_c;

    %% Store slacks
    slacks.s_diag   = s_diag;

    slacks.s_rate   = s_rate;
    slacks.s_qos    = s_qos;
    slacks.s_sig    = s_sig;
    slacks.s_int    = s_int;

    slacks.s_eh     = s_eh;
    slacks.s_rate_c = s_rate_c;
    slacks.s_sig_c  = s_sig_c;
    slacks.s_int_c  = s_int_c;
    slacks.s_qos_c  = s_qos_c;

    %% Slack report
    slack_report.total_violation = obj_prev;

    slack_report.max_diag   = max(s_diag(:));

    slack_report.max_rate   = max(s_rate(:));
    slack_report.max_qos    = max(s_qos(:));
    slack_report.max_sig    = max(s_sig(:));
    slack_report.max_int    = max(s_int(:));

    slack_report.max_eh     = max(s_eh(:));
    slack_report.max_rate_c = max(s_rate_c(:));
    slack_report.max_sig_c  = max(s_sig_c(:));
    slack_report.max_int_c  = max(s_int_c(:));
    slack_report.max_qos_c  = max(s_qos_c(:));

    %% Display summary
    disp('========== PBF FEASIBILITY SLACK REPORT ==========');
    disp(['Status              : ', status]);
    disp(['Total violation     : ', num2str(obj_prev)]);

    disp(['max s_diag          : ', num2str(slack_report.max_diag)]);
    disp(['max s_rate          : ', num2str(slack_report.max_rate)]);
    disp(['max s_qos           : ', num2str(slack_report.max_qos)]);
    disp(['max s_sig           : ', num2str(slack_report.max_sig)]);
    disp(['max s_int           : ', num2str(slack_report.max_int)]);

    disp(['max s_eh            : ', num2str(slack_report.max_eh)]);
    disp(['max s_rate_c        : ', num2str(slack_report.max_rate_c)]);
    disp(['max s_sig_c         : ', num2str(slack_report.max_sig_c)]);
    disp(['max s_int_c         : ', num2str(slack_report.max_int_c)]);
    disp(['max s_qos_c         : ', num2str(slack_report.max_qos_c)]);
    disp('==================================================');
end