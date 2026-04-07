function J_t = design_J_t(G)
    [N, K] = size(G);
    
    g_prime = (1/K) * sum(abs(G), 2); 
    
    [~, sorted_indices] = sort(g_prime, 'ascend');
    disp(size(sorted_indices));
    
    J_t = zeros(N, N);
    for i = 1:N
        J_t(i,sorted_indices(i)) = 1;
    end
end