function J_a = analyticActuatedJacobianFromQ(q, d)
% analyticActuatedJacobianFromQ  Closed-form ∂r/∂ℓ from current q
%   q = [phi1;θ1;s1; phi2;θ2;s2; phi3;θ3;s3]  (9×1)
%   d = lateral offset between chambers
%
%   Returns J_a = 3×3 Jacobian dr/dℓ.

  % Number of segments
  n = numel(q)/3;

  % Preallocate J_rq (3×3n) and J_ellq (3×3n)
  J_rq   = zeros(3,3*n);
  J_ellq = zeros(3,3*n);

  for i = 1:n
    % extract per‐segment variables
    idx    = (i-1)*3 + (1:3);
    phi    = q(idx(1));
    theta  = q(idx(2));
    s      = q(idx(3));

    %--- 1) ∂r_i/∂(φ,θ,s)  (3×3 block)  ---
    if abs(theta) < 1e-8
      % straight‐line limit → tip at [0;0;s]
      % ∂r/∂φ = [0;0;0],  ∂r/∂θ = [0;0;0],  ∂r/∂s = [0;0;1]
      Jr_block = [ 0, 0, 0;
                   0, 0, 0;
                   0, 0, 1 ];
    else
      A = 1 - cos(theta);
      B = sin(theta);
      T1 = A/theta;           % (1−cosθ)/θ
      T2 = B/theta;           % sinθ/θ
      % ∂T1/∂θ, ∂T2/∂θ
      dT1 = ( B*theta - A )/(theta^2);
      dT2 = ( cos(theta)*theta - B )/(theta^2);

      % r_i = [ s*cosφ*T1;
      %         s*sinφ*T1;
      %         s*T2 ];
      % Build partials:
      dr_dphi = [ -s*sin(phi)*T1;
                   s*cos(phi)*T1;
                   0 ];
      dr_dtheta = [  s*cos(phi)*dT1;
                     s*sin(phi)*dT1;
                     s*dT2 ];
      dr_ds = [ cos(phi)*T1;
                sin(phi)*T1;
                T2 ];

      Jr_block = [dr_dphi, dr_dtheta, dr_ds];
    end

    % insert into the big Jacobian
    J_rq(:, idx) = Jr_block;

    %--- 2) ∂ℓ/∂(φ,θ,s)  (3×3n block)  ---
    % ℓ1 = S - d*sum_k θ_k*cos(π/2 - φ_k)
    % ℓ2 = S + d*sum_k θ_k*cos(π/6 - φ_k)
    % ℓ3 = S - d*sum_k θ_k*cos(π/6 + φ_k)
    % where S = sum_k s_k

    % derivative of S wrt s_i = 1
    % ∂ℓ1/∂φ_i = -d * θ_i * sin(π/2 - φ_i) * (+1)
    % ∂ℓ1/∂θ_i = -d * cos(π/2 - φ_i)
    % ∂ℓ1/∂s_i = 1

    % similarly for ℓ2,ℓ3
    dphi1 = -d * theta * sin(pi/2 - phi);
    dtheta1 = -d * cos(pi/2 - phi);
    ds1 = 1;

    dphi2 = +d * theta * sin(pi/6 - phi);
    dtheta2 = +d * cos(pi/6 - phi);
    ds2 = 1;

    dphi3 = -d * theta * sin(pi/6 + phi);
    dtheta3 = -d * cos(pi/6 + phi);
    ds3 = 1;

    % assemble the 3×3 block for this segment
    Jell_block = [ dphi1,  dtheta1,  ds1;
                   dphi2,  dtheta2,  ds2;
                   dphi3,  dtheta3,  ds3 ];

    J_ellq(:, idx) = Jell_block;
  end

  %--- 3) chain rule:  J_a = J_rq * pinv(J_ellq)  (3×3)  ---
  % use Tikhonov regularization if you like:
  lambda = 1e-8;
  M = J_ellq.'*J_ellq + lambda*eye(3*n);
  J_a = J_rq * (M \ J_ellq.');

end