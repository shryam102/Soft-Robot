function Jxi = compute_augmented_jacobian(dh_aug)
n = size(dh_aug, 1);
T = eye(4);
p = zeros(3, n+1);
z = repmat([0;0;1], 1,n+1);

for k = 1:n
    theta = dh_aug(k,1);
    d     = dh_aug(k,2);
    a     = dh_aug(k,3);
    alpha = dh_aug(k,4);
    T = T*DH_transform(theta,d,a,alpha);
    p(:, k+1) = T(1:3,4);
    z(:, k+1) = T(1:3,3);
end
p_tip = p(:,end);
Jxi = zeros(3,n);
for k = 1:n
    idx = mod(k-1,5) + 1;
    if idx == 3
        Jxi(:,k) = z(1:3,k);
    else
        v = cross(z(:,k), p_tip - p(:,k));
        Jxi(:,k) = v(1:3);
    end
end
end