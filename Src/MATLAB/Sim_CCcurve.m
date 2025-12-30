close all; clc; clear
clear Draw_curve;
% Start and Goal Location
x_start = [10, 0, 5];
x_goal  = [5, 5, 15];

nb_steps = 2000; %maximum number of steps

x_eval = x_start;

x_path = zeros(size(x_eval));

epsilon = 1/40; %Tunable parameter
% U_path = zeros([1,1]);

for steps = 1:nb_steps
    x_path(steps,:) = x_eval;
    [U, control] = Attractive(x_goal,x_eval,2);
%     U_path(steps,:) = U;
    
    if norm(control) < 5e-2
        break;
    
    end

    x_eval = x_eval - epsilon*control;

end

for i = 1:length(x_path)
    X  = [x_path(i,1); x_path(i,2); x_path(i,3)];
%   X = [x(i);y(i); z(i)];
    IK = Inverse(X);
    pause(0.1)
    Draw_curve(IK(1), IK(2), IK(3), IK(4), x_start, x_goal);
end

if norm(x_eval - x_goal) < 5e-2
    disp('Goal Reached!');

end
