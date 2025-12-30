%% Analytic dr/dℓ derivation for 3-segment constant curvature
% Assumptions: each segment has (φ,θ,s) given by no-load map from (ℓ1,ℓ2,ℓ3)

syms ell1 ell2 ell3 d real
n = 3;              % number of identical segments

% 1) intermediate sums
S = (ell1 + ell2 + ell3)/3;

% 2) apply your no-load formulas (same φ,θ,s per segment)
phi   = atan( sqrt(3)*(ell2+ell3-2*ell1) / (3*(ell2-ell3)) );
kappa = 2*sqrt(ell1^2+ell2^2+ell3^2 - ell1*ell2 - ell2*ell3 - ell1*ell3) ...
         / (d*(ell1+ell2+ell3));
theta = S * kappa;          % total curvature
s     = S/n;                % arc‐length per segment

% 3) per‐segment tip offset r_i(φ,θ,s)
%    and sum over 3 identical segments → total r
syms i
ri = [ s*cos(phi)*(1 - cos(theta))/theta;
       s*sin(phi)*(1 - cos(theta))/theta;
            s*sin(theta)/theta ];

r = 3*ri;  % three identical contributions

% 4) compute Jacobian dr/dℓ
J_sym = simplify( jacobian(r, [ell1,ell2,ell3]) );

% 5) (optional) factor or rewrite for readability
J_simp = simplify(J_sym, 'Steps',50);

% Display result
disp('dr/dℓ ='); pretty(J_simp)
