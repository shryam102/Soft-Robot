function Draw_curve(k, theta, phi, s, Xo,Xg)
persistent endpoint_path;

T = linspace(0, theta, 100);

if k == 0
    Matrix = [linspace(0,0,100); linspace(0,0,100); linspace(0,s,100)];
else 
    r = 1/k;
    Matrix = [r*cos(phi)*(1 - cos(T)); r*sin(phi)*(1 - cos(T));  r*sin(T)];
end
plot3(Matrix(1,:), Matrix(2,:), Matrix(3,:), 'b', Linewidth = 2);
hold on
scatter3(Matrix(1,end), Matrix(2, end), Matrix(3, end), 20, 'r', 'filled');
h_start = scatter3(Xo(1), Xo(2), Xo(3), 40, [0.4660, 0.6740, 0.1880], 'filled', 'o');
h_goal =  scatter3(Xg(1), Xg(2), Xg(3), 40, [0.9290 0.6940 0.1250], 'filled', 'o');
endpoint = Matrix(:, end);
endpoint_path = [endpoint_path, endpoint];
h_path = plot3(endpoint_path(1,:), endpoint_path(2,:), endpoint_path(3,:), '--k', 'LineWidth', 1.5);

xlabel("X");
ylabel("Y");
zlabel("Z");
xlim([-20 20]);
ylim([-20,20]);
zlim([0,30]);
title("Constant Curvature Curve");
legend([h_path, h_start, h_goal], {'Path', 'Start', 'Goal'}, 'Location', 'best');
grid
drawnow, hold off





