function [li, dl_theta, dl_S] = Link_func(qi)

theta_i = qi(2);
S_i     = qi(3);

if qi(2) < 1e-3
    li = S_i;
    dl_theta = 0;
    dl_S = 1;
else
    li = 2*S_i*sin(theta_i/2)/theta_i;
    dl_theta = 2*S_i*(cos(theta_i/2)/(2*theta_i) - sin(theta_i/2)/(theta_i^2));
    dl_S = 2*sin(theta_i/2)/theta_i;

end
end