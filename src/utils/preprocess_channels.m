function channel_data = preprocess_channels(para, H_local, g_local, f_local)
    % Preprocess and organize channel data
    K = para.K; 
    K_c = para.K_c;
    scal = para.scall;
    
    % Initialize structures
    channel_data.g = cell(K, K_c);
    channel_data.g_b = cell(K, 1);
    channel_data.f = cell(K, K_c);
    channel_data.H_all = H_local * scal;
    channel_data.alpha = zeros(K, K_c);
    
    % Store channels
    for k = 1:K
        for i = 1:K_c
            channel_data.g{k,i} = g_local(:,k,i) * scal;   % N×1
            channel_data.f{k,i} = f_local(k,i);
        end
        channel_data.g_b{k} = g_local(:,k,end) * scal;  
        % channel_data.f{k,i} = f_local(k,2);     % backscatter channel
    end
    
    % Power allocation
    for k = 1:K
        for i = 1:K_c
            if i == 1
                channel_data.alpha(k,i) = para.alpha_k_n;
            else
                channel_data.alpha(k,i) = para.alpha_k_f;
            end
        end
    end
end