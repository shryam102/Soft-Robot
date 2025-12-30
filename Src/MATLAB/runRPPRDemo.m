%% runRPPRDemo.m
% Demo script to visualize inv-kinematics loop

clear; clf;
kin = RPPR_Kinematics();

% generate circle path
th=linspace(0,2*pi,100);
r=15;
x= r*cos(th); y=r*sin(th);

U = [25;30];
K = diag([.2,.2,.2,15,15,15]);

omega=zeros(3,1);
beta_prev=zeros(3,1);
L_prev=zeros(3,1);
time=0; dt=0.05;

figure; subplot(1,2,1); hold on; axis equal; grid on;
plot(x,y,'b-'); hTip=plot(0,0,'ro');
subplot(1,2,2); hold on; hL1=plot(1,U(1),'bo-'); hL2=plot(2,U(2),'go-');
axis([0 3 0 50]); grid on;

for k=1:200
    % FK no-load
    q0 = kin.q_no_load(U);
    F  = [0.2; -0.2] .* (1 - sin(0.5*time*2*pi));
    [Lnew,Bnew] = kin.compute_q_dynamics(K,q0(1:3),q0(4:6),F,omega);
    % pick next goal
    g = [x(mod(k-1,100)+1); y(mod(k-1,100)+1)];
    tip=[sum(Lnew.*sin(Bnew)); sum(Lnew.*cos(Bnew))];
    err=g-tip; J=kin.get_jacobian_actuated(U,Lnew,Bnew);
    dU = pinv(J'*J+1e-6*eye(2))*J'*err*0.2;
    U = max(min(U+dU,50),0);
    omega=(Bnew-beta_prev)/dt;

    set(hTip,'XData',tip(1),'YData',tip(2));
    set(hL1,'YData',U(1)); set(hL2,'YData',U(2));
    drawnow;
    beta_prev=Bnew; time=time+dt;
end
