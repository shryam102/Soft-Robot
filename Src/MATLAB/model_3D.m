function dq = model_3D(t,q,k,d,f_ext_fun,curves, q_ref)
n= numel(curves);
Array = cell(1,2*n);
for i=1:n
    ki = [k(i)*1.5 0; 0 k(i)];
%     delta_theta = q(2*i) - q_ref(2*i);
    di = [d(i)  0; 0 d(i)];
    Array{2*i-1} = ki;
    Array{2*i} = di;
end
K = blkdiag(Array{1:2:2*n});
D = blkdiag(Array{2:2:2*n});
J = computeSoftJacobian_3D(curves,q);
f_ext = f_ext_fun(t);
% % ---- display ----
%     disp('--------------------------------');
%     disp('  t = ');
%     disp(t);
%     disp('  J = ');
%     disp(J);
%     disp('  q = ');
%     disp(q);
%     disp('  f = ');
%     disp(f_ext);
%     disp('--------------------------------');
%     % -----------------
tau_spring = K*(q - q_ref);
    % external torque
tau_ext    = J.'*f_ext;

    % net torque into damping dynamics
tau = tau_ext - tau_spring;
dq  = D \ tau;

end