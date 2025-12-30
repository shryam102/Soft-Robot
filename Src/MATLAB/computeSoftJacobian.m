function J = computeSoftJacobian(curves,q)

L = [curves.L];

dh_aug = build_augmented_DH(curves, q);

Jxi = computePlanarAugmentedJacobian(dh_aug);

[~, Jm] = build_xi_Jm(L,q);

J = Jxi*Jm;

end