function [J2, p] = tipJacobian(L, q)
  % L: [L1 L2 L3], q: [q1; q2; q3]
  T = eye(4);
  p_all = zeros(3,4);  % p_all(:,i) = origin of frame i in base
  z_all = repmat([0;0;1],1,4);
  for i = 1:3
    Ti = trotz(q(i)) * transl(L(i),0,0);
    T  = T*Ti;
    p_all(:,i+1) = T(1:3,4);
    z_all(:,i+1) = T(1:3,3);
  end
  p_tip = p_all(:,4);
  % Build full 3×3 Jacobian
  J3 = zeros(3,3);
  for j=1:3
    J3(:,j) = cross( z_all(:,j), p_tip - p_all(:,j) );
  end
  % Return only the 2×3 translational block
  J2 = J3(1:2,:);
  p  = p_tip(1:2);
end
