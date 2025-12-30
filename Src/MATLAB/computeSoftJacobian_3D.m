function J = computeSoftJacobian_3D(curves,q)

L = [curves.L];

dh_aug = build_augmented_DH_3D(curves, q);

Jxi = computeAugmentedJacobian(dh_aug);

[~, Jm] = build_xi_Jm_3D(L,q);

J = Jxi*Jm;

end