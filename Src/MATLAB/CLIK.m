clear; close all; clc

% Initializing the Ells
u_ells = [49;45;47]; 

u_c = sum(u_ells)/50;

K = diag([15,15,.5*u_c,15,15,.5*u_c,15,15,.5*u_c]); 
D = diag([.05,.05,0.5,.05,.05,0.5,.05,.05,.5]);

goal_threshold = 0.5;

% pathCoord =  [ 0         1.8426    2.9082    2.7473    1.4278   -0.4938   -2.2072   -2.9898   -2.5115   -0.9741    0.9741    2.5115    2.9898    2.2072    0.4938   -1.4278   -2.7473   -2.9082   -1.8426  0.0;
%               -2.5000   -1.7521   -1.7357   -2.0784   -1.8660   -1.0219   -0.4228   -0.5756   -0.8526   -0.4387    0.4387    0.8526    0.5756    0.4228    1.0219    1.8660    2.0784    1.7357    1.7521  2.5;
%                8.0000    8.6494    9.2284    9.6743    9.9388    9.9932    9.8315    9.4714    8.9519    8.3292    7.6708    7.0481    6.5286    6.1685    6.0068    6.0612    6.3257    6.7716    7.3506  8.0];


% pathCoord = [0.01, -3.75, -7.5, -15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15;
%                0,    0,    0,   0,   0,   0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0,  0,  0,  0;
%               10,  7.74,  3.38,  5.48,  9.27, 11.58, 13.19, 14.35, 15.17, 15.68, 15.94, 15.94, 15.68, 15.17, 14.35, 13.19,11.58, 9.27, 5.48];

% pathCoord = [0.01, 3   3  8  11 6   8    7.5 ;
%                0,  3  -3  6  6  4  3.03   0;
%               10, 5.48 6 10  15 10 9.8019  3.38];

pathCoord = [ 8.00    7.3068    5.3644    2.5576   -0.5576   -3.3644   -5.3068   -6.0000   -5.3068   -3.3644   -0.5576    2.5576    5.3644    7.3068    8.0000
              0.00    3.0372    5.4728    6.8245    6.8245    5.4728    3.0372    0.0000   -3.0372   -5.4728   -6.8245   -6.8245   -5.4728   -3.0372   -0.0000
           10.0000    9.8019    9.2470    8.4450    7.5550    6.7530    6.1981    6.0000    6.1981    6.7530    7.5550    8.4450    9.2470    9.8019   10.0000
];

% pathCoord = [8.0;
%               12;
%               9.3];


%pointer for the goal point
numGoals = size(pathCoord,2);
i_goal = 1;
goal_coord      = pathCoord(:,i_goal);

q0_init = q_no_load(u_ells, 5, 4);
% prev_phi = q0_init(1);
% q0 = q_no_load_3(u_ells,4,4,goal_coord);
% fprintf('q0 = [%s]\n', sprintf(' %.4g', q0));
[~,~, tip_prev] = compute_Soft_curve(q0_init);


fig = figure('Name', 'IK Simulation', 'NumberTitle','off');

ax1 = subplot(1,2,1);
hold(ax1, 'on'); 
% ax1.NextPlot = 'replace';
rotate3d(ax1,'on');
grid(ax1,'on');
view(ax1,-60,50);
axis(ax1,'manual');
set(ax1, ...
      'CameraPositionMode','manual', ...
      'CameraTargetMode','manual', ...
      'CameraUpVectorMode','manual', ...
      'CameraViewAngleMode','manual');
lims = 20;
xlim(ax1,[-lims lims]);
ylim(ax1,[-lims lims]);
zlim(ax1,[0   lims]);
xlabel(ax1,'X'); ylabel(ax1,'Y'); zlabel(ax1,'Z');
title(ax1,'Backbone & Trajectory');

% initialize plot handles on ax1
path_line  = plot3(ax1, NaN,NaN,NaN,'-o','LineWidth',1);
goal_mark  = plot3(ax1, NaN,NaN,NaN,'sr','MarkerFaceColor','r');
hSkel      = plot3(ax1, NaN,NaN,NaN,'--o','LineWidth',1.5,'Color',[0 0 0 0.3]);
hSoft      = plot3(ax1, NaN,NaN,NaN,'-','LineWidth',2,'Color',[1 .5 0]);
hForce = quiver3(ax1, 0,0,0, 0,0,0, ...
                 'Color','g', ...
                 'LineWidth', 2, ...  
                 'ShowArrowHead', 'on', ...% thicker shaft
                 'Autoscale','on', ...      % larger head
                 'AutoScaleFactor',1, ...
                 'MaxHeadSize',4);
% time‐stamp text in upper‐left
hTime = text(0.02,0.95,'Time: 0.00 s', ...
                 'Units','normalized','FontSize',12, ...
                 'FontWeight','bold','Parent',ax1);

ax2 = subplot(1,2,2);
hold(ax2,'on');
grid(ax2,'on');
xlim(ax2,[0 4]); ylim(ax2,[0 max(u_ells)*1.1]);
xlabel(ax2,'Chamber #'); ylabel(ax2,'Length (mm)');
title(ax2,'u\_ells');
len1_mark = plot(ax2,[1 1],[0 u_ells(1)],'-o','LineWidth',3);
len2_mark = plot(ax2,[2 2],[0 u_ells(2)],'-o','LineWidth',3);
len3_mark = plot(ax2,[3 3],[0 u_ells(3)],'-o','LineWidth',3);

drawnow;

% time_prev = tic;
% t_elapsed = 0;

dt = 0.02;
t_elapsed = 0;


% f_ext  = @(t) [0; 0; 0.3];
% f_ext =  @(t) [0; 0; -.3*sin(2*pi*0.2*t)];
% f_ext =  @(t) [0; 0; -1*abs(sin(2*pi*0.2*t))];
f_ext =  @(t) [0; 0; -1];

Kp     = 1.5;     
gain_min = 0.2;  
gain_max = 3.0;   
tol_eq   = 0.1;  




while ishandle(fig)
%     tNow = toc(time_prev);
%     dt = tNow;
%     time_prev = tic;
    goal_coord      = pathCoord(:,i_goal);
    
    tStart = tic;

    t_elapsed = t_elapsed+ dt;
    
    hTime.String = sprintf('Time: %.2f s', t_elapsed);


    q0 = q_no_load(u_ells,3,4);

%     prev_phi = q0(1);
%     q0 = q_no_load_3(u_ells,4,4,goal_coord);
    disp('-----------------------------------------------');
    fprintf('U_ells: [%s]\n', sprintf(' %.4g', u_ells));
    fprintf('Time_elapsed: [%s]\n', sprintf(' %.4g', t_elapsed ));
    fprintf('dt = [%s]\n', sprintf(' %.4g', dt));
    fprintf('q0 = [%s]\n', sprintf(' %.4g', q0));
    q_now = q_dynamics(q0, f_ext, t_elapsed, dt, K, D);
    fprintf('qnow = [%s]\n', sprintf(' %.4g', q_now));
    fprintf('dt = [%s]\n', sprintf(' %.4g', dt));

    [softPts,skelPts,tip_coord] = compute_Soft_curve(q_now);
    fvec = f_ext(t_elapsed);
    error = goal_coord - tip_coord;
    tip_vel = (tip_coord - tip_prev)/dt;
    tip_prev = tip_coord;
    fprintf('Goal_coord = [%s]\n', sprintf(' %.4g', goal_coord));
    fprintf('Tip_coord  = [%s]\n', sprintf(' %.4g', tip_coord));

    % ------ CLIK ---------------

    J_ac = Actuated_Jacobian(q_now);
    [U,S,V] = svd(J_ac, 'econ');
    sigma = diag(S);
    sigma_0 = 0.01*max(sigma);
    nu = 50;
    h = (sigma.^3 + nu*sigma.^2 + 2*sigma + 2*sigma_0) ...
    ./ (sigma.^2 + nu*sigma + 2);
%     lambda = 1e-6;
%     M = J_ac.'*J_ac + lambda*eye(size(J_ac,2));
%     J_ac_inv = M\J_ac.';
    J_ac_inv = V * diag(1./h) * U';
    rate_fb = 0.09 * (1 + max(0, 2 - norm(error)));
    rate_ff = 5e-6;
    delta_u = J_ac_inv*(rate_fb*error + rate_ff*(-tip_vel));
    fprintf('delta_u = [%s]\n', sprintf(' %.4g', delta_u));

%     gain = compute_gain(u_ells, error, Kp, gain_min, gain_max, tol_eq);
    gain = 0.5;
    u_ells = u_ells + gain*delta_u;
    
%     u_ells = min(max(u_ells + real(delta_u), 0), 50);
    u_ells = min(max(u_ells,0),50);
%     tolerance = 0.05; % Define what "close" means
%     if (abs(u_ells(1) - u_ells(2)) < tolerance && ...
%         abs(u_ells(2) - u_ells(3)) < tolerance && ...
%         abs(u_ells(1) - u_ells(3)) < tolerance)
%     
%     % Add small random noise to each element
%         noise_magnitude = 0.1; % Adjust as needed
%         noise = (rand(3,1) - 0.5) * 2 * noise_magnitude; % Random noise between -0.05 and 0.05
%         u_ells = u_ells + noise;
%     
%     % Ensure bounds are still respected after adding noise
%         u_ells = min(max(u_ells,0),50);
%     end

    u_ells = round(u_ells,3);

    fprintf('u_ells = [%s]\n', sprintf(' %.4g', u_ells));
    set(len1_mark, 'YData', [0,u_ells(1)]);
    set(len2_mark, 'YData', [0,u_ells(2)]);
    set(len3_mark, 'YData', [0,u_ells(3)]);

    if norm(error) < goal_threshold
        i_goal = mod(i_goal,numGoals) + 1;
    end

    set(path_line, 'XData',pathCoord(1,:), 'YData', pathCoord(2,:), 'ZData', pathCoord(3,:));
    set(goal_mark, 'XData',goal_coord(1,:),'YData', goal_coord(2,:),'ZData', goal_coord(3,:));
    set(hSkel, 'XData', skelPts(1,:), 'YData', skelPts(2,:), 'ZData', skelPts(3,:));
    set(hSoft, 'XData', softPts(1,:), 'YData', softPts(2,:), 'ZData', softPts(3,:));
    set(hForce, ...
        'XData', tip_coord(1), 'YData', tip_coord(2), 'ZData', tip_coord(3), ...
        'UData', fvec(1)*2, ...
        'VData', fvec(2)*2, ...
        'WData', fvec(3)*2);

   drawnow;
   pause(0.02);
%     compTime = toc(tStart);
%     pauseTime = max(0, dt - compTime);
%     pause(pauseTime);

end

function keyPressCallback(src, event)
    switch event.Key
        case 'space'
            % Toggle pause state
            current_state = getappdata(src, 'pause_state');
            setappdata(src, 'pause_state', ~current_state);
            
            if ~current_state
                fprintf('Simulation PAUSED - Press SPACEBAR to resume\n');
            else
                fprintf('Simulation RESUMED\n');
            end
            
        case 'escape'
            % Close the figure to exit simulation
            fprintf('Simulation terminated by user\n');
            close(src);
    end
end


function [li, dl_theta, dl_S] = Link_func(qi)

theta_i = qi(2);
S_i     = qi(3);

if qi(2) < 1e-3
    li = S_i;
    dl_theta = 0;
    dl_S = 1;
else
    li = 2*S_i*sin(theta_i/2)/theta_i;
    dl_theta = 2*S_i*(cos(theta_i/2)/(2*theta_i) - sin(theta_i/2)/(theta_i^2));
    dl_S = 2*sin(theta_i/2)/theta_i;

end
end

function [xi, Jmi] = single_map(qi)
phi   = qi(1);
theta = qi(2);

[li, dl_theta, dl_S] = Link_func(qi);

xi = [   phi;
      theta/2;
          li;
      theta/2;
        -phi];

Jmi = [1      0       0;
       0     1/2      0;
       0   dl_theta  dl_S;
       0     1/2      0;
      -1      0       0];

end

function [xi, Jm] = build_XI_J(q)
n = numel(q)/3;
xi = zeros(5*n,1);
Jm = zeros(5*n,3*n);
for i = 1:n
    idx = (i-1)*5 + (1:5);
    [xi_block, Jmi] = single_map(q(3*i-2:3*i));
    xi(idx) = xi_block;
    Jm(idx, 3*i-2:3*i) = Jmi;
end

end

function dh_aug = augmented_dh(q)

n = numel(q)/3;
dh_aug = zeros(5*n, 4);
[Xi, ~] = build_XI_J(q);

for i = 1:n
    xi = Xi((i-1)*5 + (1:5));
    base = (i - 1)*5;
    dh_aug(base +1, :) = [xi(1),   0,    0, -pi/2];
    dh_aug(base +2, :) = [xi(2),   0,    0,  pi/2];
    dh_aug(base +3, :) = [0      xi(3),  0  -pi/2];
    dh_aug(base +4, :) = [xi(4),   0,    0,  pi/2];
    dh_aug(base +5, :) = [xi(5),   0,    0,     0];
end
end

function Jxi = compute_augmented_jacobian(dh_aug)
n = size(dh_aug, 1);
T = eye(4);
p = zeros(3, n+1);
z = repmat([0;0;1], 1,n+1);

for k = 1:n
    theta = dh_aug(k,1);
    d     = dh_aug(k,2);
    a     = dh_aug(k,3);
    alpha = dh_aug(k,4);
    T = T*DH_transform(theta,d,a,alpha);
    p(:, k+1) = T(1:3,4);
    z(:, k+1) = T(1:3,3);
end
p_tip = p(:,end);
Jxi = zeros(3,n);
for k = 1:n
    idx = mod(k-1,5) + 1;
    if idx == 3
        Jxi(:,k) = z(1:3,k);
    else
        v = cross(z(:,k), p_tip - p(:,k));
        Jxi(:,k) = v(1:3);
    end
end
end

function J = computeJacobian(q)
dh_aug = augmented_dh(q);
Jxi = compute_augmented_jacobian(dh_aug);

[~, Jm] = build_XI_J(q);

J = Jxi*Jm;

end

function q_next = q_dynamics(q_0,f_ext,t0,dt, K, D)
tspan = t0: dt: t0 + dt;

[~,Q] = ode45(@(t,q) model_3(t,q,f_ext,K,D, q_0), tspan, q_0); 
q_next = Q(end,:)';
end

function [softPts, skelPts, r] = compute_Soft_curve(q)
n = numel(q)/3;
samples = 100;
T = eye(4);
softPts = zeros(3,samples*n);
skelPts = zeros(3, n+1);
phis = q(1:3:end);
thetas = q(2:3:end);
L = q(3:3:end);

for i =1:n
    Xi = position_vector(thetas(i), samples, phis(i), L(i), T);
    idx = (i-1)*samples + (1:samples);
    softPts(:,idx) = Xi;
    T = T*Transform(thetas(i), phis(i), L(i));
    skelPts(:, i+1) = softPts(:,samples*i);

end
r = softPts(:,end);
end

function J_a = Actuated_Jacobian(q)
% n = numel(q)/3;
% lambda = 1e-6;
J_rq  = computeJacobian(q);
J_Eq  = compute_Jacobian_Ell_Q(q);
J_Eq_inv = pinv(J_Eq);

% M = J_Eq.'* J_Eq + lambda*eye(3*n);
% J_Eq_reg = M\J_Eq.';

J_a = J_rq * J_Eq_inv;

end

function g = compute_gain(ell, err, Kp, g_min, g_max, tol_eq)
   
    if abs(ell(1)-ell(2))<tol_eq && abs(ell(2)-ell(3))<tol_eq
        g = g_min;
    else
        e = norm(err);
        g = Kp * e;
        g = min(max(g, g_min), g_max);
    end
end
