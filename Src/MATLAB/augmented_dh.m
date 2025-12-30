function dh_aug = augmented_dh(q)

n = numel(q)/3;
dh_aug = zeros(5*n, 4);
[Xi, ~] = build_XI_J(q);

for i = 1:n
    xi = Xi((i-1)*5 + (1:5));
    base = (i - 1)*5;
    dh_aug(base +1, :) = [xi(1),   0,    0, -pi/2];
    dh_aug(base +2, :) = [xi(2),   0,    0,  pi/2];
    dh_aug(base +3, :) = [0      xi(3),  0  -pi/2];
    dh_aug(base +4, :) = [xi(4),   0,    0,  pi/2];
    dh_aug(base +5, :) = [xi(5),   0,    0,     0];
end
end