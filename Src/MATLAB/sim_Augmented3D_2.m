%% sim_Augmented3D.m
clear, clc, close all

%--- setup curves and parameters ------------------------------------------
curves    = curve_params();  
n         = numel(curves);
tspan     = 0:.1:30;
q0        = repmat([0; pi/3],n,1);      % [0;pi/3;0;pi/3;…]
k_val     = 5;
d_val     = 10;

%--- precompute the big K matrix once ------------------------------------
% each 2×2 Ki = [0 0; 0 k_val]
Kcells = repmat({ [0 0; 0 k_val] }, 1, n);
Kbig   = blkdiag( Kcells{:} );

%--- external‐force handle ------------------------------------------------
f_ext_fun_3D = @(t) [0.4*sin(2*pi*0.2*t); 0; 0];

%--- ODE options & solve --------------------------------------------------
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t1, q1] = ode45( ...
    @(t,q) fastModel3D(t,q,Kbig,d_val,f_ext_fun_3D,curves), ...
    tspan, q0, opts );

% … plot / animate q1 as before …
