function FW_Curve(theta, phi, l, color)

n = length(theta);
T = eye(4);
samples = 100;
X_all = zeros(3, samples*n); 

for i = 1:n
    theta_i = theta(i);
    phi_i = phi(i);
    l_i = l(i);
    Xi = position_vector(theta_i,samples,phi_i,l_i,T);
    idx = (i-1)*samples + 1 : i*samples;  
    X_all(:, idx) = Xi;
    T = T*Transform(theta_i,phi_i,l_i);

end

plot3(X_all(1,:), X_all(2,:), X_all(3,:), 'color', color, Linewidth = 2);

