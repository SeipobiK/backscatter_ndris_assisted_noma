function Theta = initialize_ris_phases(para)
    % Initialize RIS phase shift matrix with random phases
    u = exp(1i * pi * (2 * rand(para.N, 1)));
    Theta = diag(u);
end