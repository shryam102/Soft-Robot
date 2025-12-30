function q_0 = q_no_load(ell, n, d)

q_0 = zeros(3*n,1);

ell_1 = ell(1);
ell_2 = ell(2); 
ell_3 = ell(3);

S = (ell_1 + ell_2 + ell_3)/3;

 phi = atan2(sqrt(3) * (ell_2 + ell_3 - 2*ell_1), (3*(ell_2 - ell_3)));
 kappa = 2*sqrt(ell_1^2 + ell_2^2 + ell_3^2 - ell_1 * ell_2 - ell_1*ell_3 - ell_2*ell_3)/(d*(ell_1 + ell_2 + ell_3));
 theta = S * kappa;




if (abs(ell_1-ell_2) <= 0.001) && (abs(ell_2-ell_3)<=0.001)
    phi = 0;
    kappa = 1/1000000000;
    theta = S * kappa;
    
% else
%    
%     if phi < -0.01
%     phi = 2*pi - abs(phi);
%     end
end

% if (abs(ell_2 + ell_3 - 2*ell_1) < 1e-4) && (ell_2 > ell_3)
%     phi = 0;
% elseif (abs(ell_2 + ell_3 - 2*ell_1) < 1e-4) && (ell_2 < ell_3)
%     phi = pi;
% 
% elseif (abs(ell_2 - ell_3) < 1e-4) && (ell_2 + ell_3 > 2*ell_1)
%     phi = pi/2;
% 
% elseif (abs(ell_2 - ell_3) < 1e-4) && (ell_2 + ell_3 < 2*ell_1)
%     phi = 3*pi/2;
% 
% 
% 
% end
% 
% if (phi < 0) && (ell_2 < ell_3)
%     phi = pi - abs(phi);
% elseif (phi > 0) && (ell_2 < ell_3)
%     phi = pi + abs(phi);
% elseif (phi < 0) && (ell_2 + ell_3 < 2*ell_1)
%     phi = 2*pi - abs(phi);
% 
% end




for i = 1:3:3*n
    q_0(i) = phi;
    q_0(1 + i) = theta/n;
    q_0(2 + i) = S/n;
end
end