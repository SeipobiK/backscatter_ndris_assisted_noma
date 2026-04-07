function [V_max,max_eigenvalue]=max_eigVect(V_opt)
        [V,D] = eig(V_opt);
        [max_eigenvalue, index] = max(diag(D));
        V_max =V(:, index); 
end