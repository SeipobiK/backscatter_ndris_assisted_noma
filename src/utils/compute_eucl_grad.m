function [grad] = compute_eucl_grad(para,w_k,Theta,channel_data,decoding_order,...
    beta,zeta,y,z,duals_vars,duals_vars_c,alpha,eta)

%% Parameters
K     = para.K;
K_c   = para.K_c;

gamma_noma = 2^para.R_min_n - 1;
gamma_back = 2^para.R_min_c - 1;

%% Initialize Euclidean gradient (NxN)
grad = zeros(size(Theta));

%% MAIN LOOP
for k = 1:K

    order_k = decoding_order(k,:);
    strong_user = order_k(end);

    for i = 1:K_c

        %% =========================
        % Desired cascaded channel
        %% =========================
        d_ki = channel_data.g{k,i}' * Theta * channel_data.H_all * w_k(:,k);

        %% =========================
        % Gradient of B(k,i)
        %% =========================
        grad_B = zeros(size(Theta));

        % Inter-cluster interference
        for j = 1:K
            if j ~= k
                d_inter = channel_data.g{k,i}' * Theta * channel_data.H_all * w_k(:,j);
                grad_B = grad_B + 2*d_inter * conj(channel_data.g{k,i}) * (channel_data.H_all * w_k(:,j))';
            end
        end

        % Intra cluster interference
        pos = find(order_k == i);
        for idx = pos+1:K_c
            j_user = order_k(idx);
            grad_B = grad_B + 2*alpha(k,j_user)*d_ki * conj(channel_data.g{k,i}) * (channel_data.H_all * w_k(:,k))';
        end

        % Backscatter interference (FIXED: use current cluster k)
        for j = 1:K
            d_b = channel_data.g_b{k}' * Theta * channel_data.H_all * w_k(:,j);
            grad_B = grad_B + 2*eta(j)*d_b * conj(channel_data.g_b{k}) * (channel_data.f{k,i} * channel_data.H_all * w_k(:,j))';
        end

        %% =========================
        % QT NOMA gradient
        %% =========================
        grad = grad + ...
            y(k,i)*sqrt((1+beta(k,i))*alpha(k,i)) * conj(channel_data.g{k,i}) * (channel_data.H_all * w_k(:,k))' ...
            - abs(y(k,i))^2 * ( ...
                2*alpha(k,i)*d_ki * conj(channel_data.g{k,i}) * (channel_data.H_all * w_k(:,k))' ...
                + grad_B );

        %% =========================
        % NOMA QoS dual gradient
        %% =========================
        grad = grad + duals_vars(k,i) * ( ...
            2*alpha(k,i)*d_ki * conj(channel_data.g{k,i}) * (channel_data.H_all * w_k(:,k))' ...
            - gamma_noma * grad_B );

        %% =========================
        % Backscatter decoder (strong user only)
        %% =========================
        if i == strong_user

            d_b = channel_data.g_b{k}' * Theta * channel_data.H_all * channel_data.f{k,i} * w_k(:,k);

            %% Bc gradient
            grad_Bc = zeros(size(Theta));
            for j = 1:K
                if j ~= k
                    d_bc = channel_data.g_b{k}' * Theta * channel_data.H_all * channel_data.f{k,i} * w_k(:,j);
                    grad_Bc = grad_Bc + 2*d_bc * conj(channel_data.g_b{k}) * (channel_data.f{k,i} * channel_data.H_all * w_k(:,j))';
                end
            end

            %% QT backscatter term (FIXED indexing: z(k), zeta(k))
            grad = grad + ...
                z(k)*sqrt((1+zeta(k))*eta(k)) * conj(channel_data.g_b{k}) * (channel_data.f{k,i} * channel_data.H_all * w_k(:,k))' ...
                - abs(z(k))^2 * ( ...
                    2*eta(k)*d_b * conj(channel_data.g_b{k}) * (channel_data.f{k,i} * channel_data.H_all * w_k(:,k))' ...
                    + grad_Bc );

            %% Backscatter QoS dual (FIXED indexing)
            grad = grad + duals_vars_c(k) * ( ...
                2*eta(k)*d_b * conj(channel_data.g_b{k}) * (channel_data.f{k,i} * channel_data.H_all * w_k(:,k))' ...
                - gamma_back * grad_Bc );
        end
    end
end

% %% Symmetry penalty gradient
% if isfield(para, 'nu')
%     nu = para.nu;
%     grad = grad - 4*nu*(Theta - Theta.');
% end

end