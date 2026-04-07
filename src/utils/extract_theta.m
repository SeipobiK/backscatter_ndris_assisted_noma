function [Theta_opt, bestWSR] = extract_theta( ...
    V_opt, para, H, H_c,decoding_order,alpha,w_k)

    % Step 1: Dominant eigenvector extraction
    [V_max, max_eigenvalue_v] = max_eigVect(V_opt);
    v_k = sqrt(max_eigenvalue_v) * V_max;

    % Step 2: Generate candidate phase vectors
    cand = { ...
        exp(1j * angle(v_k)), ...
        exp(-1j * angle(v_k)), ...
        conj(exp(1j * angle(v_k))), ...
        conj(exp(-1j * angle(v_k))) ...
    };
    Theta_opt = diag(cand{2}); % Default to first candidate (can be changed based on evaluation)
    % Step 3: Evaluate candidates
    % bestWSR = -inf;
    % Theta_opt = [];

    % for t = 1:length(cand)
    %     theta_test = diag(cand{t});

    %     [WSR,~,~,~,~,~,~,~,~,~,~] = compute_wsr(...
    %     para, H, H_c,decoding_order,alpha,w_k);

    %     if WSR > bestWSR
    %         bestWSR = WSR;
    %         Theta_opt = theta_test;
    %         % disp(['Best t: ', num2str(t)]);
    %     end
    % end
end