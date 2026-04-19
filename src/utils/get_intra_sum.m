function intra_sum = get_intra_sum(alpha, decoding_order, K, K_c)
    intra_sum = zeros(K, K_c);
    for k = 1:K
        for i = 1:K_c
            sum_alpha = 0;
            for j = 1:K_c
                if decoding_order(k, j) > decoding_order(k, i)
                    sum_alpha = sum_alpha + alpha(k, decoding_order(k, j));
                end
            end
            intra_sum(k,i) = sum_alpha;
        end
    end
end