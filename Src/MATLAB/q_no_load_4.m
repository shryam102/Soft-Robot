function q_0 = q_no_load_4(ell, n, d, phi_prev)
    q_0 = zeros(3*n,1);
    ell_1 = ell(1);
    ell_2 = ell(2);
    ell_3 = ell(3);
    S = (ell_1 + ell_2 + ell_3)/3;
    
    % Calculate phi as before
    phi = atan2(sqrt(3) * (ell_2 + ell_3 - 2*ell_1), (3*(ell_2 - ell_3)));
    kappa = 2*sqrt(ell_1^2 + ell_2^2 + ell_3^2 - ell_1 * ell_2 - ell_1*ell_3 - ell_2*ell_3)/(d*(ell_1 + ell_2 + ell_3));
    theta = S * kappa;
    
    if (abs(ell_1-ell_2) <= 0.001) && (abs(ell_2-ell_3)<=0.001)
        phi = 0;
        kappa = 1/1000000000;
        theta = S * kappa;
    end
    
    % Apply phi continuity constraint if phi_prev is provided
    if nargin >= 4 && ~isempty(phi_prev)
        max_phi_change = deg2rad(10); % 20 degrees in radians
        
        % Calculate the difference, considering periodic nature of angles
        phi_diff = phi - phi_prev;
        
        % Wrap the difference to [-pi, pi]
        phi_diff = atan2(sin(phi_diff), cos(phi_diff));
        
        % If the difference exceeds the maximum allowed change, constrain it
        if abs(phi_diff) > max_phi_change
            if phi_diff > 0
                phi = phi_prev + max_phi_change;
            else
                phi = phi_prev - max_phi_change;
            end
            
            % Wrap phi to [0, 2*pi] or [-pi, pi] as needed
            phi = atan2(sin(phi), cos(phi));
            if phi < 0
                phi = phi + 2*pi;
            end
        end
    end
    
    % Populate q_0 as before
    for i = 1:3:3*n
        q_0(i) = phi;
        q_0(1 + i) = theta/n;
        q_0(2 + i) = S/n;
    end
end