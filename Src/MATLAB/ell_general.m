function ell = ell_general(q)
phi   = q(1:3:end);
theta = q(2:3:end);
s     = q(3:3:end);
d = 4;
Ssum = sum(s);

ell_1 = Ssum - d*sum(theta .*sin(phi));
ell_2 = Ssum + d*sum(theta .*cos(pi/6 - phi));
ell_3 = Ssum - d*sum(theta .*cos(pi/6 + phi));
ell = [ell_1;ell_2;ell_3];
end