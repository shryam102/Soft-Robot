function [n, dn_t, dn_l] = N_theta(L, q2)
    [B, dB_t, dB_l] = B_theta(L, q2);

 
    x = (L/(q2 * B)) * sin(q2/2);

    x = max(-1, min(1, x));

   
    num = 0.5*cos(q2/2)*(q2*B) -   sin(q2/2)*(B + q2*dB_t);
    dx_t = L * num / ((q2 * B)^2);
    dx_l = (sin(q2/2)/(q2*B)) - L*sin(q2/2)*dB_l/(q2*B^2); 


    if abs(q2) < 1e-3
        n   = 0;
        dn_t = 0;
        dn_l = 0;
        return;
    end

    n = acos(x);

    tol = 1e-8;
    if abs(1 - abs(x)) < tol
        dn_t = 0;
        dn_l = 0;
    else
        dn_t = -dx_t / sqrt(1 - x^2);
        dn_l = -dx_l / sqrt(1 - x^2);
    end
end
