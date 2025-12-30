function q_next = q_dynamics(q_0,f_ext,t0,dt, K, D)
tspan = t0: dt: t0 + dt;

[~,Q] = ode45(@(t,q) model_3(t,q,f_ext,K,D, q_0), tspan, q_0); 
q_next = Q(end,:)';
end
