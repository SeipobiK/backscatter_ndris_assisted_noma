function [BS_array, RIS_array] = generate_arrays(para)

    % Author: Xidong Mu (base) | Modified by Kgomotjo Seipobi
    BS_array = [(1:para.M)'-(para.M+1)/2, zeros(para.M,1), zeros(para.M,1)];
    
    RIS_array = zeros(para.N, 3);
    for i = 1:para.RIS_size(1)
        for  j = 1:para.RIS_size(2)
            n = (i-1)*para.RIS_size(2) + j;
            RIS_array(n,:) = [ i, 0, j ];
        end
    end
    RIS_array = RIS_array- [(para.RIS_size(1)+1)/2, 0, 0];
    
end