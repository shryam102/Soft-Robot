function Draw_curves(Cell, curves)
n = numel(Cell);

L = curves.L;

draw_multi_segment(Cell{1}, L, 'r');
hold on;
draw_multi_segment(Cell{2}, L, 'b');
FW_Curve(Cell{3}, 0, L*n, 'g');
xlabel("X");
ylabel("Y");
zlabel("Z");
xlim([-L*n,L*n]);
ylim([-L*n,L*n]);
zlim([0,L*n]);
grid on;
drawnow; 
hold off
