function dq = model(t, q,k,d,f_ext_fun,curves)
n = numel(curves);
k_vals = linspace(k,k,n);
d_vals = linspace(d,d,n);

K = diag(k_vals(1:n));
D = diag(d_vals(1:n));

J = -computeSoftJacobian(curves,q);
f_ext = f_ext_fun(t);
tau = J.' * f_ext;
dq = D \ (tau - K*q);
end