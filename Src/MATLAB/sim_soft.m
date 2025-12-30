clear all; close all; clc

%Material Properties
m = 1;
d = 5;
d_theta = 5;
K = 30;
Ks = 10;
l = 1;

total_time = 15; %Time for simulation

F = -10;  % Constant X direction Force
J = -10; % Impulse in X Direction 

t_impulse = 3; % Time at which impulse is provided

y0 = [l; 0; 0; 0];  % Initial State

tspan1 = 0:0.1:t_impulse;  
tspan2 = t_impulse:0.1:total_time; 

[t1, y1] = ode45(@(t,y) soft(y, m, d, K, Ks, F,l, d_theta), tspan1, y0);

y0_impulse = y1(end, :); 
y0_impulse(2) = y0_impulse(2) - J*sin(y0_impulse(3))/ m; % Apply impulse to velocity
y0_impulse(4) = y0_impulse(4) - J*cos(y0_impulse(3))/ (m * (m*y0_impulse(1)^2));

[t2, y2] = ode45(@(t,y) soft(y, m, d, K, Ks, F,l, d_theta), tspan2, y0_impulse);


t = [t1; t2];
y = [y1; y2];

for k=1:length(t)
    pause(0.03)
    draw_soft(y(k,:));
end

figure;

subplot(2,1,1);
plot(t, y(:,1),'red',Linewidth = 1.5);
title('Variation of Length');
xlabel('Time');
ylabel('Length')


subplot(2,1,2);
plot(t, y(:,3), 'blue', LineWidth=1.5);
title("Variation of Theta vs Time");
xlabel('Time');
ylabel('Theta');

grid

