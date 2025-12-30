function dq = model_cc(t,q,K,D,f_ext_fun,curves)
J = -jacobian_cc(q,curves);
f_ext = f_ext_fun(t);
tau = J' * f_ext;

dq = (1/D)*(tau - K*q);

end

