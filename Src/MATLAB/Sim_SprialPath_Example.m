clc; clear; close all;

clear Draw_curve

z = linspace(5,25,200);

t = linspace(0,4*pi,200);

x = 15*cos(t);
y = 15*sin(t);


for i = 1:length(t)
%     X  = [path(i,1); path(i,2); path(i,3)];
    X = [x(i);y(i); z(i)];
    IK = Inverse(X);
    pause(0.01)
    Draw_curve(IK(1), IK(2), IK(3), IK(4), [x(1), y(1), z(1)], [x(end), y(end), z(end)]);
end

