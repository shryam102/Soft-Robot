%% fastModel3D.m
function dq = fastModel3D(t, q, Kbig, d, f_ext_fun, curves)
  % Number of 2×2 segments:
  n = numel(curves);

  % 1) compute generalized torque
  J   = computeSoftJacobian_3D(curves, q);  % your Jacobian
  tau = J' * f_ext_fun(t);

  % 2) compute RHS = tau – K*q
  rhs = tau - Kbig*q;

  % 3) allocate output
  dq = zeros(size(q));

  % 4) fill in dq block by block (element-wise division)
  %    block i has di = [d1 0; 0 d2], so
  %      dq(2i-1) = rhs(2i-1)/d1
  %      dq(2i)   = rhs(2i)  /d2
  %    where d1 = d*(q(2i)^2) + 1      (your baseline reg = 1)
  for i = 1:n
    idx1 = 2*i-1;
    idx2 = idx1 + 1;

    d1 = d * ( q(idx2)^2 ) + 1;  % baseline eps = 1
    d2 = d;

    dq(idx1) = rhs(idx1) / d1;
    dq(idx2) = rhs(idx2) / d2;
  end
end
