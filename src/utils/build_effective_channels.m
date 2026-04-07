function [H_eff, Hc_eff] = build_effective_channels(para, channel_data, Theta,J_r,J_t)
    % Build effective channels for all users and subcarriers
    K = para.K; 
    K_c = para.K_c;
    
    H_eff = cell(K, K_c);
    Hc_eff = cell(K, K_c);
    
    % Design J matrices (assuming these functions exist)
    % You may need to adjust these based on your implementation
    % g_LOS_reshaped = reshape(channel_data.g{1,1}, para.N, []);
    % J_r = design_J_r(g_LOS_reshaped);
    % J_t = design_J_t(channel_data.H_all);
    
    for k = 1:K
        for i = 1:K_c
            % Direct link: H = g' * J_t' * Theta' * J_r' * H_all
            % disp(channel_data.g{k,i});
            % disp(Theta);
            H_eff{k,i} = channel_data.g{k,i}' * J_r * Theta * J_t * channel_data.H_all;
            
            % Backscatter link: H_c = g_b' * J_t' * Theta' * J_r' * H_all * f
            Hc_eff{k,i} = channel_data.g_b{k}' * J_r * Theta * J_t * ...
                          channel_data.H_all * channel_data.f{k,i};
        end
    end
end