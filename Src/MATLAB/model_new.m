function dq = model_new(t, q, K, D, f_ext_fun, curves)
    J     = -1*computeSoftJacobian(curves, q);
    f_ext = f_ext_fun(t);
    
    % ---- display ----
    disp('--------------------------------');
    disp('  t = ');
    disp(t);
    disp('  J = ');
    disp(J);
    disp('  q = ');
    disp(q);
    disp('  f = ');
    disp(f_ext);
    disp('--------------------------------');
    % -----------------
    
    tau = J.' * f_ext;
    dq  = D \ (tau - K*q);
end
