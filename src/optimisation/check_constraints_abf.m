function [W_opt, A_opt, B_opt, A_c_opt, B_c_opt, obj_prev, status, slack_report, slacks, duals] = ...
    check_constraints_abf(para, channel_data, H, H_c, ...
    A_prev, B_prev, A_c_prev, B_c_prev, decoding_order, alpha, eta)

    %% Parameters
    K       = para.K;
    K_c     = para.K_c;
    M       = para.M;
    P_max   = para.P_max;
    noise   = para.noise;
    R_min   = para.R_min_n;
    R_min_c = para.R_min_c;
    eh      = para.bst_threshold;
    rho     = para.rho;

    %% CVX optimization
    cvx_begin quiet
        cvx_solver mosek

        %% Variables
        variable W(M,M,K) hermitian semidefinite

        variable A(K,K_c) nonnegative
        variable B(K,K_c) nonnegative
        variable A_c(K) nonnegative
        variable B_c(K) nonnegative

        variable R(K,K_c) nonnegative
        variable R_c(K) nonnegative

        %% Slack variables
        variable s_pow nonnegative

        variable s_rate(K,K_c) nonnegative
        variable s_qos(K,K_c) nonnegative
        variable s_sig(K,K_c) nonnegative
        variable s_int(K,K_c) nonnegative

        variable s_eh(K) nonnegative
        variable s_rate_c(K) nonnegative
        variable s_sig_c(K) nonnegative
        variable s_int_c(K) nonnegative
        variable s_qos_c(K) nonnegative

        %% Dual variables
        dual variable d_pow

        dual variables d_rate{K,K_c}
        dual variables d_qos{K,K_c}
        dual variables d_sig{K,K_c}
        dual variables d_int{K,K_c}

        dual variables d_eh{K}
        dual variables d_rate_c{K}
        dual variables d_sig_c{K}
        dual variables d_int_c{K}
        dual variables d_qos_c{K}

        %% Objective: minimise total violation
        minimize( ...
            s_pow + ...
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

            %% Power constraint
            total_power = 0;

            for tx_cluster = 1:K
                total_power = total_power + real(trace(W(:,:,tx_cluster)));
            end

            d_pow : total_power <= P_max + s_pow;

            %% Main constraints
            for user_cluster = 1:K

                order_k = decoding_order(user_cluster,:);
                strong_user = order_k(end);

                for user_idx = 1:K_c

                    %% NOMA SCA rate approximation
                    rate_sca = log2(1 + 1/(A_prev(user_cluster,user_idx) * B_prev(user_cluster,user_idx))) ...
                        - (log2(exp(1)) / ...
                        (A_prev(user_cluster,user_idx) * ...
                        (1 + A_prev(user_cluster,user_idx) * B_prev(user_cluster,user_idx)))) ...
                        * (A(user_cluster,user_idx) - A_prev(user_cluster,user_idx)) ...
                        - (log2(exp(1)) / ...
                        (B_prev(user_cluster,user_idx) * ...
                        (1 + A_prev(user_cluster,user_idx) * B_prev(user_cluster,user_idx)))) ...
                        * (B(user_cluster,user_idx) - B_prev(user_cluster,user_idx));

                    d_rate{user_cluster,user_idx} : ...
                        R(user_cluster,user_idx) <= ...
                        rate_sca + s_rate(user_cluster,user_idx);

                    %% NOMA QoS
                    d_qos{user_cluster,user_idx} : ...
                        R(user_cluster,user_idx) + ...
                        s_qos(user_cluster,user_idx) >= R_min;

                    %% Inter-cluster interference
                    inter = 0;
                    inter_b = 0;

                    for interferer_cluster = 1:K
                        if interferer_cluster ~= user_cluster

                            H_user = H{user_cluster,user_idx}' * ...
                                     H{user_cluster,user_idx};

                            Hc_user = H_c{user_cluster,user_idx}' * ...
                                      H_c{user_cluster,user_idx};

                            inter = inter + real(trace( ...
                                W(:,:,interferer_cluster) * H_user));

                            inter_b = inter_b + eta(interferer_cluster) * ...
                                real(trace(W(:,:,interferer_cluster) * Hc_user));
                        end
                    end

                    %% Intra-cluster interference
                    intra = 0;

                    pos_user = find(order_k == user_idx);
                    H_user = H{user_cluster,user_idx}' * H{user_cluster,user_idx};

                    for order_pos = pos_user+1:K_c
                        interfering_user = order_k(order_pos);

                        intra = intra + alpha(user_cluster,interfering_user) * ...
                            real(trace(W(:,:,user_cluster) * H_user));
                    end

                    %% NOMA signal constraint
                    signal_noma = alpha(user_cluster,user_idx) * ...
                        real(trace(W(:,:,user_cluster) * H_user));

                    d_sig{user_cluster,user_idx} : ...
                        inv_pos(A(user_cluster,user_idx)) <= ...
                        signal_noma + s_sig(user_cluster,user_idx);

                    %% NOMA interference constraint
                    Hc_user = H_c{user_cluster,user_idx}' * ...
                              H_c{user_cluster,user_idx};

                    self_backscatter_int = eta(user_cluster) * ...
                        real(trace(W(:,:,user_cluster) * Hc_user));

                    d_int{user_cluster,user_idx} : ...
                        B(user_cluster,user_idx) + ...
                        s_int(user_cluster,user_idx) >= ...
                        intra + inter + inter_b + ...
                        self_backscatter_int + noise;

                    %% Backscatter constraints only for strong/SIC user
                    if user_idx == strong_user

                        %% EH channel
                        % H_c{k,i} = h_b{k,i} * f{k,i}
                        % therefore h_b{k,i} = H_c{k,i}/f{k,i}
                        h_eh = H_c{user_cluster,user_idx} / ...
                               channel_data.f{user_cluster,user_idx};

                        H_eh = h_eh' * h_eh;

                        harvested_energy = (1 - eta(user_cluster)) * rho * ...
                            real(trace(W(:,:,user_cluster) * H_eh));

                        d_eh{user_cluster} : ...
                            harvested_energy + s_eh(user_cluster) >= eh;

                        %% Backscatter signal constraint
                        signal_c = eta(user_cluster) * ...
                            real(trace(W(:,:,user_cluster) * Hc_user));

                        d_sig_c{user_cluster} : ...
                            inv_pos(A_c(user_cluster)) <= ...
                            signal_c + s_sig_c(user_cluster);

                        %% Backscatter SCA rate approximation
                        rate_c_sca = log2(1 + 1/(A_c_prev(user_cluster) * B_c_prev(user_cluster))) ...
                            - (log2(exp(1)) / ...
                            (A_c_prev(user_cluster) * ...
                            (1 + A_c_prev(user_cluster) * B_c_prev(user_cluster)))) ...
                            * (A_c(user_cluster) - A_c_prev(user_cluster)) ...
                            - (log2(exp(1)) / ...
                            (B_c_prev(user_cluster) * ...
                            (1 + A_c_prev(user_cluster) * B_c_prev(user_cluster)))) ...
                            * (B_c(user_cluster) - B_c_prev(user_cluster));

                        d_rate_c{user_cluster} : ...
                            R_c(user_cluster) <= ...
                            rate_c_sca + s_rate_c(user_cluster);

                        %% Backscatter interference constraint
                        d_int_c{user_cluster} : ...
                            B_c(user_cluster) + s_int_c(user_cluster) >= ...
                            inter + inter_b + noise;

                        %% Backscatter QoS
                        d_qos_c{user_cluster} : ...
                            R_c(user_cluster) + ...
                            s_qos_c(user_cluster) >= R_min_c;
                    end
                end
            end

    cvx_end

    %% Outputs
    obj_prev = cvx_optval;
    status   = cvx_status;

    W_opt   = W;
    A_opt   = A;
    B_opt   = B;
    A_c_opt = A_c;
    B_c_opt = B_c;

    %% Slacks
    slacks.s_pow    = s_pow;

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

    slack_report.max_pow    = max(s_pow(:));

    slack_report.max_rate   = max(s_rate(:));
    slack_report.max_qos    = max(s_qos(:));
    slack_report.max_sig    = max(s_sig(:));
    slack_report.max_int    = max(s_int(:));

    slack_report.max_eh     = max(s_eh(:));
    slack_report.max_rate_c = max(s_rate_c(:));
    slack_report.max_sig_c  = max(s_sig_c(:));
    slack_report.max_int_c  = max(s_int_c(:));
    slack_report.max_qos_c  = max(s_qos_c(:));

    %% Duals
    duals.d_pow    = d_pow;

    duals.d_rate   = d_rate;
    duals.d_qos    = d_qos;
    duals.d_sig    = d_sig;
    duals.d_int    = d_int;

    duals.d_eh     = d_eh;
    duals.d_rate_c = d_rate_c;
    duals.d_sig_c  = d_sig_c;
    duals.d_int_c  = d_int_c;
    duals.d_qos_c  = d_qos_c;

    %% Display summary
    disp('========== ABF FEASIBILITY SLACK REPORT ==========');
    disp(['Status              : ', status]);
    disp(['Total violation     : ', num2str(obj_prev)]);

    disp(['max s_pow           : ', num2str(slack_report.max_pow)]);
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