function J = Jacobian_ell_q(q)
n = numel(q)/3;
h = 1e-6;
e0 = ell_general(q);
J = zeros(3,3*n);
for i = 1:3*n
    dq = zeros(3*n, 1);
    dq(i) = h;
    J(:,i) = (ell_general(q + dq) - e0)/h;
end
end