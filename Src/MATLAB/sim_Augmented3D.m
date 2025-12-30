clear all, close all,  clc
curves = curve_params();
n = numel(curves);

tspan = 0:.05:30;
q0_1 = repmat([0; 0],n,1);
% q0_2 = zeros(n,1);
% q0_cc = 0;
k_val = 5;
d_val = 10;



% f_ext_fun = @(t) [0.4*sin(2*pi*0.2*t); 0];
f_ext_fun_3D = @(t) [-1.4*sin(2*pi*0.2*t); 0; 0];
% f_ext_fun = @(t) [1 + t/100; 0];
% f_ext = [1;0];
opts = odeset('RelTol',1e-6, 'AbsTol',1e-6);
[t_1,q_1] = ode45(@(t,q) model_3D(t,q,k_val,d_val,f_ext_fun_3D,curves), tspan,q0_1);
% [t_2,q_2] = ode45(@(t,q) model(t,q,k_val,d_val,f_ext_fun,curves), tspan,q0_2);
% [t_3,q_3] = ode45(@(t,q) model_cc(t,q,k_val,d_val,f_ext_fun,curves), tspan,q0_cc);

% len = size(q_2,1);
% q_2 = [zeros(len,1), q_2(:,1), zeros(len,1), q_2(:,2), zeros(len,1), q_2(:,3)];





% then continue the animation
for k = 2:length(t_1)
    draw_multi_segment(q_1(k,:)', curves, 'r');
%     pause(0.05);
end




