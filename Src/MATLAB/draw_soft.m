function draw_soft(state)

x = state(1);
th = state(3);

pendx = -x*sin(th);
pendy = x*cos(th);


plot([-10, 10], [0,0], 'k', Linewidth = 1), hold on
plot([0 pendx],[0 pendy],'k','LineWidth',2); 
rectangle('Position', [-0.1, -0.1, 2*0.1, 2*0.1], 'Curvature' , 1, 'FaceColor' ,[1 .1 .1], 'LineWidth', .5);
axis([-1 1 -2 2]);
axis equal
drawnow, hold off


