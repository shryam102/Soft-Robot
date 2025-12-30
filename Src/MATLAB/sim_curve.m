clear all, close all, clc
curves = curve_params();
n = numel(curves);

tspan = 0:.1:30;
q0_1 = zeros(n,1);
q0_cc = 0;
k_val = 5;
d_val = 10;

k_vals = [5,5,5];
d_vals = [10,10,10];

K = diag(k_vals(1:n));
D = diag(d_vals(1:n));

f_ext_fun = @(t) [0.4*sin(2*pi*0.2*t); 0];
% f_ext_fun = @(t) [1 + t/100; 0];
% f_ext = [1;0];
opts = odeset('RelTol',1e-6, 'AbsTol',1e-6);
[t_1,q_1] = ode45(@(t,q) model(t,q,K,D,f_ext_fun,curves), tspan,q0_1,opts);
[t_2,q_2] = ode45(@(t,q) model_cc(t,q,k_val,d_val,f_ext_fun,curves), tspan,q0_cc,opts);

phi = [0;0;0];
phi_cc = 0;

% for k=1:length(t_1)
%     q_list = q_1(k,:)';
%     Q = [0;q_list(1); 0;q_list(2);0;q_list(2)];
% %     Q_cc = [phi_cc,q_2(k)];
%     L = curves.L;
%     draw_multi_segment(Q,L);
%     pause(0.05);
% end

