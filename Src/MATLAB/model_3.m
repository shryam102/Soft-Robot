function dq = model_3(t,q,f_ext_fun, K,D, q_ref)
n= numel(q)/3;
% Array = cell(1,2*n);
% for i=1:n
%     ki = [k(i)*1.5 0; 0 k(i)];
% %     delta_theta = q(2*i) - q_ref(2*i);
%     di = [d(i)  0; 0 d(i)];
%     Array{2*i-1} = ki;
%     Array{2*i} = di;
% end
% K = diag([15,10,.2,10,7,.2,10,7,.2]);
% D = diag([.05,.05,0.05,.05,.05,0.05,.05,.05,.05]);
J = computeJacobian(q);
f_ext = f_ext_fun(t);
% % ---- display ----
%     disp('--------------------------------');
%     disp('  t = ');
%     disp(t_elapsed);
% %     disp('  J = ');
% %     disp(J);
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