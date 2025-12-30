function dx = soft(x, m, d, K, Ks, F,l,d_theta)


Sx = sin(x(3));

Cx = cos(x(3));

I = m*x(1)^2;

%System of Equations

dx(1,1) = x(2);
dx(2,1) = (1/m)*(-F*Sx - K*(x(1) - l) - d*x(2) + m*x(1)*x(4)^2);
dx(3,1) = x(4);
dx(4,1) = (1/I)*(-3*F*x(1)*Cx - 3*Ks*x(3) - 2*m*x(1)*x(2)*x(4) - d_theta*x(4));

