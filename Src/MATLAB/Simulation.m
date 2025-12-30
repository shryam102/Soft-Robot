clear; close all; clc
curves = curve_params();
n      = 4;
tspan  = 0:.05:30;
K = diag([15,10,.2,10,7,.2,10,7,.2,10,7,.2]); 
D = diag([.05,.05,0.05,.05,.05,0.05,.05,.05,.05,.05,.05,.05]);
q0     = repmat([pi/100;pi/9; 8.33],n,1);
f_ext  = @(t) (t<=20 && t>=5)*[0; 0; -0.4*sin(2*pi*0.2*t)];
% f_ext = [0.2;0;-0.2];
opts   = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t_1, q_1] = ode45(@(t,q) model_3(t,q,f_ext,K,D,q0), tspan, q0);

%% One-Time axes & Camera setup-

figure;
ax = gca;
ax.NextPlot = 'replacechildren';
rotate3d(ax,'on');
grid(ax,'on');
view(ax, -60, 50);
axis(ax,'manual');
set(ax, ...
  'CameraPositionMode',  'manual', ...
  'CameraTargetMode',    'manual', ...
  'CameraUpVectorMode',  'manual', ...
  'CameraViewAngleMode', 'manual');

L    = 10;
lims = L * n;
xlim(ax,[-lims-5, lims+5]);
ylim(ax,[-(lims+5), lims+5]);
zlim(ax,[0,      lims]);
xlabel('X'); ylabel('Y'); zlabel('Z');
drawnow;

%% Pre-create line, arrow & time text

hSkel  = plot3(ax, NaN,NaN,NaN, '--o','LineWidth',1.5,'Color',[0 0 0 0.3]);
hold(ax,'on');
hSoft  = plot3(ax, NaN,NaN,NaN, '-',  'LineWidth',2,  'Color','r');
hForce = quiver3(ax, 0,0,0, 0,0,0, ...
                 'Color','b', ...
                 'LineWidth', 2, ...  
                 'ShowArrowHead', 'on', ...% thicker shaft
                 'Autoscale','on', ...      % larger head
                 'AutoScaleFactor',1, ...
                 'MaxHeadSize',4);
% time display in upper left
hTime  = text(0.02, 0.95, 'Time: 0.00 s', ...
              'Units','normalized', ...
              'FontSize',12, ...
              'FontWeight','bold', ...
              'Color','k');
hold(ax,'off');
drawnow;

%% Animate (30 s real‐time)
startTime = tic;
scale     = L/3;   % arrow length scale
for k = 1:length(t_1)
    qk      = q_1(k,:)';
%     skelPts = computeSkeleton(qk, curves);
%     softPts = computeSoftCurve(qk, curves);

    [softPts, skelPts, ~] = compute_Soft_curve(qk);

    % update skeleton & curve
    set(hSkel, 'XData', skelPts(1,:), 'YData', skelPts(2,:), 'ZData', skelPts(3,:));
    set(hSoft, 'XData', softPts(1,:), 'YData', softPts(2,:), 'ZData', softPts(3,:));

    % update force arrow at the tip
    tip  = skelPts(:,end);
    fvec = f_ext(t_1(k));
    set(hForce, ...
        'XData', tip(1), 'YData', tip(2), 'ZData', tip(3), ...
        'UData', fvec(1)*5, ...
        'VData', fvec(2)*5, ...
        'WData', fvec(3)*5);

    % update time text
    hTime.String = sprintf('Time: %.2f s', t_1(k));

    drawnow limitrate

    % pace to 30 s total
    elapsed = toc(startTime);
    pause(max(t_1(k) - elapsed, 0));
end
