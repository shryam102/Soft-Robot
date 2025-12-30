function J = jacobian_cc(q,curves)

L = [curves.L];
L_cc = length(L)*L(1);

if abs(q) < 1e-8
   half_chord = L_cc/2;
else
   half_chord = (L_cc*sin(q/2))/q;
end

dh_pars     = [-q/2,    0,         0,      -pi/2; 
                0  half_chord      0,        0;
                0, half_chord      0,       pi/2;
              -q/2,     0,         0,         0];

n = size(dh_pars,1);

T = eye(4);

for i = 1:n
    theta = dh_pars(i,1);
    d     = dh_pars(i,2);
    a     = dh_pars(i,3);
    alpha = dh_pars(i,4);
    T = T*DH_transform(theta, d, a, alpha);
end

J = cross(T(1:3,3), T(1:3, 4));
J = J(1:2,1);
