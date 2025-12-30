function dh_aug = build_augmented_DH(curves, q)
n = numel(curves);
dh_aug = zeros(4*n, 5);

for i = 1:n
    L   = curves(i).L;
    mu  = curves(i).mu;
    qi  = q(i);
    base = (i-1)*4;
    if abs(qi) < 1e-8
        half_chord = L/2;
    else
        half_chord = (L*sin(qi/2))/qi;
    end
    dh_aug(base+1, :) = [ -qi/2,      0,         0,   -pi/2,  0 ];
    dh_aug(base+2, :) = [ 0,      half_chord,    0,      0,  mu ];
    dh_aug(base+3, :) = [ 0,      half_chord,    0,    pi/2,  0 ];
    dh_aug(base+4, :) = [ -qi/2,      0,         0,      0,  0 ];
end
end
