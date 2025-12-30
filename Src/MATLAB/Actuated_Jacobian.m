function J_a = Actuated_Jacobian(q)
n = numel(q)/3;
lambda = 1e-6;
J_rq  = computeJacobian(q);
J_Eq  = compute_Jacobian_Ell_Q(q);

M = J_Eq.'* J_Eq + lambda*eye(3*n);
J_Eq_reg = M\J_Eq.';

J_a = J_rq * J_Eq_reg;

end
