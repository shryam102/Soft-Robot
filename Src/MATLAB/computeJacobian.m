function J = computeJacobian(q)
dh_aug = augmented_dh(q);
Jxi = compute_augmented_jacobian(dh_aug);

[~, Jm] = build_XI_J(q);

J = Jxi*Jm;

end