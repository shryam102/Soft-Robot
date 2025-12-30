function dh_aug = build_augmented_DH_3D(curves, q)
n = numel(curves);
dh_aug = zeros(10*n, 5);
[Xi,~] = build_xi_Jm_3D([curves.L],q);

for i = 1:n
    mu  = curves(i).mu;
    xi = Xi((i-1)*10 + (1:10));
    base = (i-1)*10;
    dh_aug(base+1, :) = [ xi(1),      0,         0,   -pi/2,  0 ];
    dh_aug(base+2, :) = [ xi(2),      0,         0,    pi/2,  0 ];
    dh_aug(base+3, :) = [   0,      xi(3),       0,   -pi/2,  0 ];
    dh_aug(base+4, :) = [ xi(4),      0,         0,    pi/2,  0 ];
    dh_aug(base+5, :) = [ xi(5),      0,         0,      0,   mu];
    dh_aug(base+6, :) = [ xi(6),      0,         0,   -pi/2,  0 ];
    dh_aug(base+7, :) = [ xi(7),      0,         0,    pi/2,  0 ];
    dh_aug(base+8, :) = [   0,      xi(8),       0,   -pi/2,  0 ];
    dh_aug(base+9, :) = [ xi(9),      0,         0,    pi/2,  0 ];
    dh_aug(base+10,:) = [ xi(10),     0,         0,      0,   0 ];  
end
end
