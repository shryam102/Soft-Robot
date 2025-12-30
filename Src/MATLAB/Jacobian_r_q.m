function J = Jacobian_r_q(q)
n = numel(q)/3;
h  = 1e-6;
[~,~,r0] = compute_Soft_curve(q);
J = zeros(3,3*n);
for i = 1:3*n
    dq = zeros(3*n,1);
    dq(i) = h;
    [~,~, r] = compute_Soft_curve(q + dq);
    J(:,i) = (r - r0)/h;
end
end