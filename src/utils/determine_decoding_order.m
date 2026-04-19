function [decoding_order,gains_it] = determine_decoding_order(para, H_eff, W)
    % Determine SIC decoding order based on channel gains
    K = para.K;
    K_c = para.K_c;
    gains_it= zeros(K, K_c); % Store effective channel gains for each user
 
    decoding_order = zeros(K, K_c);
    
    for k = 1:K
        gains = zeros(1, K_c);
        if K_c == 2

                for i = 1:K_c
                    h = H_eff{k,i};
                                inter = 0;
                                inter_b = 0;

                                for j = 1:K
                                    if j ~= k
                                        inter   = inter  + abs(h * W(:, j))^2;
                                        inter_b = inter_b + abs(h * W(:, j))^2 * para.eta;
                                    end
                                end
                    
                    gains(i) = abs(h * W(:, k))^2/(inter + para.noise);  % Effective channel gain considering interference and noise
                    gains_it(k,i)= gains(i);
                end
        elseif K_c == 4

                for i = 1:K_c
                    
                    h = H_eff{k,i};
                                inter = 0;
                                inter_b = 0;

                                for j = 1:K
                                    if j ~= k
                                        inter   = inter  + abs(h * W(:, j))^2;
                                        inter_b = inter_b + abs(h * W(:, j))^2 * para.eta;
                                    end
                                end
                    gains(i) = abs(h * W(:, k))^2/(inter + para.noise);  % Effective channel gain considering interference 
                    gains_it(k,i)= gains(i);
                end

        end
        % Sort in ascending order (weakest first)
        [~, order] = sort(gains, 'ascend');
        decoding_order(k, :) = order;
    end  
    % Display first two decoding orders for verification (optional)
    if K >= 2 && nargout == 0
        fprintf('Decoding order - User 1: %s\n', mat2str(decoding_order(1,:)));
        fprintf('Decoding order - User 2: %s\n', mat2str(decoding_order(2,:)));
    end
end