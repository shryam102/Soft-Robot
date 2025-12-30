function [x,d_x] = f_theta(theta)


if abs(theta) < 1e-3
    x = 0;
    d_x = 0;

else
    x = (theta^2 - 2*cos(theta) -2*theta*sin(theta) + 2)/(theta)^2;
    d_x = (4*theta*sin(theta) - 2*(theta^2)*cos(theta) + 4*(cos(theta) - 1))/(theta^3);

end



