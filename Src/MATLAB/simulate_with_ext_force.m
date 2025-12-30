% simulate_force_pcc.m
% Simulate and visualize a 3-segment PCC soft robot under tip force

function simulate_force_pcc()
    %--- Robot parameters ---
    L = [0.1, 0.1, 0.1];      % segment lengths (m)
    n = numel(L);

    %--- Compliance gains ---
    K = diag([50, 50, 50]);     % stiffness (Nm/rad)
    D = diag([5,  5,  5 ]);     % damping   (Nms/rad)

    %--- External tip force ---
    f_ext_fun = @(t) [ 1.0*sin(2*pi*0.5*t);   % 0.5 Hz sine in X
                      -0.5*cos(2*pi*0.2*t)]; % 0.2 Hz cosine in Y

    %--- ODE definition: D qdot + K q = J(q)' * f_ext  --> qdot = D^{-1}(J' f_ext - K q)
    odefun = @(t, q) D \ ( tipJacobian(L, q)'*f_ext_fun(t) - K*q );

    %--- Simulation ---
    tspan = [0, 5];             % simulate 5 seconds
    q0    = zeros(n,1);         % start from zero curvature
    opts  = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [T, Q] = ode45(odefun, tspan, q0, opts);

    %--- Set up figure ---
    figure('Color','w'); axis equal;
    totalL = sum(L);
    xlim([-totalL, totalL]); ylim([-totalL, totalL]);
    xlabel('X (m)'); ylabel('Y (m)');
    title('PCC Soft Robot under Tip Force'); drawnow;

%     %--- Animate ---
%     for k = 1:5:length(T)
%         qk = Q(k, :)';   % current curvatures
%         drawPCC(L, qk);
%         drawnow;
%         if k < length(T)
%             cla;         % clear axes for next frame
%         end
%     end
% end

%% Helper: planar Jacobian of tip (2×3)
function J2 = tipJacobian(L, q)
    % Build end-effector Jacobian for 3 revolute joints in plane
    T = eye(4);
    p_all = zeros(3,4);
    z_all = repmat([0;0;1], 1, 4);
    for i = 1:3
        T = T * trotz(q(i)) * transl(L(i), 0, 0);
        p_all(:,i+1) = T(1:3,4);
        z_all(:,i+1) = T(1:3,3);
    end
    p_tip = p_all(:,4);
    J3 = zeros(3,3);
    for j = 1:3
        J3(:,j) = cross(z_all(:,j), p_tip - p_all(:,j));
    end
    J2 = J3(1:2, :);
end

%% Helper: draw PCC shape in 2D
function drawPCC(L, q)
    pts = [0;0];
    theta_acc = 0;
    for i = 1:numel(L)
        Li = L(i); qi = q(i);
        if abs(qi) < 1e-3
            xi = [0, Li]';
        else
            R = Li/qi;
            u = linspace(0,1,30);
            xi = [ R*sin(qi*u); R*(1-cos(qi*u)) ];
        end
        Racc = [cos(theta_acc), -sin(theta_acc);
                sin(theta_acc),  cos(theta_acc)];
        Xi   = (Racc * xi) + pts;
        plot(Xi(1,:), Xi(2,:), 'b-', 'LineWidth', 2);
        pts = Xi(:,end);
        theta_acc = theta_acc + qi;
        hold on;
    end
    plot(pts(1), pts(2), 'ro', 'MarkerSize',8, 'MarkerFaceColor','r'); % tip
    hold off;
end

%% Basic DH helpers
function T = trotz(theta)
    T = [cos(theta), -sin(theta), 0, 0;
         sin(theta),  cos(theta), 0, 0;
               0    ,        0   , 1, 0;
               0    ,        0   , 0, 1];
end

function T = transl(x,y,z)
    T = eye(4);
    T(1:3,4) = [x; y; z];
end
