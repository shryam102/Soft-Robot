clear; clc; close all;

%% Plotting
f = figure;
set(f, 'Color', 'w', 'Position', [926, 80, 560, 420])
hold on
axis equal
grid on
title('Tip Least Squares Annealing')

%% Initialization
% param = init_sba_param();
param = [];

% model tip frame
p0 = [0, pi/9, 3]';
[~, ~, g0] = compute_Soft_curve(p0);
plotTransforms(se3(g0), FrameLabel="model")

% actual tip frame
pd = p0 + rand(size(p0));
[~, ~, gd] = compute_Soft_curve(pd);
plotTransforms(se3(gd), FrameLabel="cam")

%% Gradient descent hyperparameters
tol = 0.1; % error tolerance
max_itr = 200; % max iterations
epsilon = 1e-10; % finite difference perturbation
s0 = 1; % default step size

%% Gradient descent
itr = 0;
err = dist_func(g0, gd); 
fprintf('Initial error:\t%.2f\n', err);

p_cur = p0;
% tic
while err > tol && itr < max_itr
    tic
    a = grad(p_cur, gd, epsilon, param);
    s = armijo(s0, p_cur, a, gd, param);

    p_cur = p_cur - s*a;
    toc

    [~, ~, g_cur] = compute_Soft_curve(p_cur);
    err = dist_func(g_cur, gd);
    % err_list = [err_list, err];
    itr = itr + 1;
end
% err_list
% toc

plotTransforms(se3(g_cur), FrameLabel="corrected")
fprintf('Final error:\t%.2f\n', err);
%% Gradient descent functions

function a = grad(pc, g_des, epsilon, param)
n = 3;
X = eye(3);

a = zeros(n, 1);
[~, ~, g0] = compute_Soft_curve(pc);
d0 = dist_func(g0, g_des);

for i = 1:n
    pl = pc + epsilon*X(:, i);
%     gl = pcc_param_to_T(pl(1), pl(2), pl(3), param);
    [~, ~, gl] = compute_Soft_curve(pl);
    dl = dist_func(gl, g_des);

    %% Central difference
    % pr = pc - epsilon*X(:, i);
    % gr = pcc_param_to_T(pr(1), pr(2), pr(3), param);
    % dr = dist_func(gr, g_des);

    % a(i) = (dl - dr)/(2*epsilon);

    %% Forward difference
    a(i) = (dl - d0)/epsilon;
end
end

function d = dist_func(g1, g2)
% returns distance d between two SE3 elements based on definition
%% Pose
d = sqrt(trace((g1 - g2)*transpose(g1 - g2)));

%% Position
% d = norm(g1(1:3, 4) - g2(1:3, 4));
end

function s = armijo(s0, pc, a, gd, param)
sigma = 0.5; % 0 < sigma < 1
beta = 0.75; % 0 < beta < 1


s = s0;

% gc = pcc_param_to_T(pc(1), pc(2), pc(3), param);
[~, ~, gc] = compute_Soft_curve(pc);
fc = dist_func(gc, gd);
while true
    pn = pc - s*a;
%     gn = pcc_param_to_T(pn(1), pn(2), pn(3), param);
    [~, ~, gn] = compute_Soft_curve(pn);
    fn = dist_func(gn, gd);
    if fn <= fc - sigma*s*a'*a
        return
    else
        s = beta*s;
    end
end
end