function [a] = steering_vector(array,theta, phi)
    % Author: Xidong Mu (base) | Modified by Kgomotjo Seipobi

    a = exp(-1i*array*K(theta,phi));
    end
    
    function k = K(theta,phi)
    k = pi * [cos(theta).*cos(phi), sin(theta).*cos(phi), sin(phi)]';
end
    