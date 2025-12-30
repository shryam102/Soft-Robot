% part2.m
clear; close all; clc

%% 1) Define plant & weights
A  = [1 2;
      0 1];
B  = [0;1];
Bw = [1;1];          % disturbance enters both states equally
Cz = [1 1];          % performance output z = x1 + x2

Q = Cz'*Cz;          % = [1 1; 1 1]
R = 1;               % cost on u^2

%% 2) H2 (LQR) static state‐feedback
[P,~,~] = care(A, B, Q, R);
K_H2    = -B' * P;   % 1×2 row vector

%% 3) H∞ static state‐feedback for gamma = 6,8,10
gammas = [6, 8, 10];
K_Hinf = zeros(numel(gammas), 2);

for i = 1:numel(gammas)
    gamma = gammas(i);

    % Build generalized plant Pgen with outputs [ z ; y = x ] stacked
    %  - first row is Cz*x  (performance)
    %  - next two rows are identity*x (full-state measurement)
    Pgen = ss( A, [Bw B], ...
               [Cz; eye(2)], ...
               zeros(3,2) );

    % Synthesize a (static) H∞ controller via Riccati
    %   nmeas = 2 (we measure both x1,x2), nctrl = 1
    %   'METHOD','ric' enforces the ARE solution
    [Kdyn,~,~] = hinfsyn(Pgen, 2, 1, 'METHOD','ric', 'GMIN', gamma);

    % Extract the direct‐feedthrough matrix D ∈ ℝ^(1×2)
    % which is exactly the static gain K = -B' * X
    [~,~,~,Dk] = ssdata(Kdyn);

    K_Hinf(i,:) = Dk;
end

%% 4) Simulate w(t) = sin(2πt) on [0,10]
t    = linspace(0,10,5000);
w    = sin(2*pi*t);
x0   = [0;0];
Z2   = zeros(numel(gammas), numel(t));
Zh   = zeros(numel(gammas), numel(t));

for i = 1:numel(gammas)
    % H2 closed-loop
    sys2 = ss(A + B*K_H2, Bw, Cz, 0);
    Z2(i,:) = lsim(sys2, w, t, x0)';

    % H∞ closed-loop
    sysh = ss(A + B*K_Hinf(i,:)', Bw, Cz, 0);
    Zh(i,:) = lsim(sysh, w, t, x0)';
end

%% 5) Plot results
figure('Position',[100 100 600 500]);
for i = 1:numel(gammas)
    subplot(numel(gammas),1,i);
    plot(t, Z2(i,:), 'b', t, Zh(i,:), 'r--','LineWidth',1.2);
    title(sprintf('\\gamma = %d', gammas(i)));
    ylabel('z(t)');
    if i==numel(gammas)
        xlabel('Time (s)');
    end
    legend('H_2','H_\infty','Location','NorthWest');
    grid on;
end
