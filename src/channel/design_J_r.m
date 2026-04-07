function J_r = design_J_r(G)
    [N, M] = size(G);

    g_prime = (1/M) * sum(abs(G), 2);
    [~, sorted_indices] = sort(g_prime, 'ascend');
    
    J_r = zeros(N, N);
    for i = 1:N
        J_r(sorted_indices(i),i) = 1;
    end
end