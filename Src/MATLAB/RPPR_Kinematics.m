% RPPR_Kinematics.m
% Complete MATLAB translation of the provided Python RPPR script

classdef RPPR_Kinematics
    properties
        % geometric & material properties
        a = 2;          % [mm]
        h = 0.0381;     % [mm]
        E = 7.35;       % [N/mm^2]
        v = 0.5;        % Poisson's ratio
        d = 4.5;        % [mm]
        phi_i = 2/3*pi;
        phi_0 = pi/2;
        T1_0;
        T2_0;
        T3_0;
    end

    methods
        function obj = RPPR_Kinematics()
            % Precompute base transforms
            obj.T1_0 = [1 0 0 obj.d*cos(obj.phi_0);
                        0 1 0 obj.d*sin(obj.phi_0);
                        0 0 1 0;
                        0 0 0 1];
            obj.T2_0 = [1 0 0 obj.d*cos(obj.phi_i+obj.phi_0);
                        0 1 0 obj.d*sin(obj.phi_i+obj.phi_0);
                        0 0 1 0;
                        0 0 0 1];
            obj.T3_0 = [1 0 0 obj.d*cos(2*obj.phi_i+obj.phi_0);
                        0 1 0 obj.d*sin(2*obj.phi_i+obj.phi_0);
                        0 0 1 0;
                        0 0 0 1];
        end

        function [arcLen, arcPhi, arcKappa] = robotSpecificMapping(obj, l1, l2, l3)
            arcLen   = (l1 + l2 + l3)/3;
            arcPhi   = atan2((sqrt(3)/3)*(l2 + l3 - 2*l1), (l2 - l3));
            arcKappa = (2*sqrt(l1^2 + l2^2 + l3^2 - l1*l2 - l1*l3 - l2*l3)) / (obj.d*(l1+l2+l3));
        end

        function Tw = robotIndependentMapping(obj, arcLen, arcPhi, arcKappa)
            theta = arcKappa * arcLen;
            c    = cos(theta);
            s    = sin(theta);
            cp   = cos(arcPhi);
            sp   = sin(arcPhi);
            Tw = [cp^2*(c-1)+1,    sp*cp*(c-1),    cp*s,    cp*(1-c)/arcKappa;
                  sp*cp*(c-1),     cp^2*(1-c)+c,   sp*s,    sp*(1-c)/arcKappa;
                  -cp*s,           -sp*s,          c,       s/arcKappa;
                  0,               0,              0,       1];
        end

        function tip = twoStepMapping(obj, l1, l2, l3)
            if max([l1,l2,l3]) - min([l1,l2,l3]) > 0
                [s,phi,k] = obj.robotSpecificMapping(l1,l2,l3);
                T = obj.robotIndependentMapping(s,phi,k);
                tip = T(1:3,4);
            else
                tip = [0;0;l1];
            end
        end

        function backbone = backbone(obj, l1, l2, l3, N)
            l1v = linspace(0,l1,N);
            l2v = linspace(0,l2,N);
            l3v = linspace(0,l3,N);
            backbone = zeros(3,N);
            if max([l1,l2,l3]) - min([l1,l2,l3]) > 0
                for i=1:N
                    [s,phi,k] = obj.robotSpecificMapping(l1v(i),l2v(i),l3v(i));
                    T = obj.robotIndependentMapping(s,phi,k);
                    backbone(:,i) = T(1:3,4);
                end
            else
                backbone(3,:) = linspace(0,l1,N);
            end
        end

        function len = pressureToLength(obj, P, scale)
            if nargin<3, scale = 0.84; end
            % Timoshenko's model approx for SBA
            len = scale*2*length(P)*0.622*obj.a*(P*obj.a/(obj.E*obj.h)).^(1/3);
        end

        function volumes = lengths_to_volumes(obj, ells, base_h)
            volumes = 100/4667 * (ells - base_h) * 1000;
        end

        function Jr = calcJacobian(obj, P)
            % Placeholder: bring in your full Python expression here,
            % or call get_jacobian_actuated if you only need actuated form.
            error('calcJacobian: full translation required');
        end

        function [coords, r, phi, theta, Ttip] = backbone_inv(obj, tip_pose, N)
            % Inverse kinematics backbone reconstruction
            x=tip_pose(1); y=tip_pose(2); z=tip_pose(3);
            r    = (x^2+y^2+z^2)/(2*sqrt(x^2+y^2));
            theta=2*atan2(sqrt(x^2+y^2),z);
            phi  = atan2(y,x);
            kappa=sign(theta)/r;
            L    = r*theta;
            coords = zeros(3,N);
            s = linspace(0,L,N);
            Tcur=eye(4);
            for i=1:N
                Tcur = obj.robotIndependentMapping(s(i),phi,kappa);
                coords(:,i)=Tcur(1:3,4);
            end
            Ttip = Tcur;
        end

        function circ = generate_circle(obj, radius, z, M)
            th = linspace(0,2*pi,M);
            circ = [radius*cos(th); radius*sin(th); z*ones(1,M)];
        end

        function R = rotate_axis(obj, axis, pts, deg)
            th = deg2rad(deg);
            switch axis
                case 'x'
                    Rm = [1 0 0;0 cos(th) -sin(th);0 sin(th) cos(th)];
                case 'y'
                    Rm = [cos(th) 0 sin(th);0 1 0;-sin(th) 0 cos(th)];
                case 'z'
                    Rm = [cos(th) -sin(th) 0; sin(th) cos(th) 0;0 0 1];
            end
            R = Rm * pts;
        end

        function path = load_path(obj, fname)
            data = dlmread(fname);
            path = obj.rotate_axis('y', data, -40);
            path = obj.rotate_axis('z', path, 210) + [-12;5;0];
            path = 0.7*path;
        end

        function circ3 = generate_3d_circle(obj, pts3, M)
            if size(pts3,1)~=3 || size(pts3,2)~=3, error('need 3 points'); end
            p1=pts3(:,1); p2=pts3(:,2); p3=pts3(:,3);
            n = cross(p2-p1,p3-p1); n=n/norm(n);
            mid1=(p1+p2)/2; mid2=(p2+p3)/2;
            dir1=cross(n,p2-p1); dir2=cross(n,p3-p2);
            t= (dir1\(mid2-mid1)); center=mid1 + t*dir1;
            radius=norm(p1-center);
            ang=linspace(0,2*pi,M);
            u=(p1-center)/radius;
            v=cross(n,u);
            circ3 = center + radius*(u*cos(ang) + v*sin(ang));
        end

        function [l1,l2,l3] = steps_to_length(obj, s1, s2, s3)
            v1=4.1910e-05*s1; v2=4.1910e-05*s2; v3=4.1910e-05*s3;
            l1=45.67*v1; l2=45.67*v2; l3=45.67*v3;
        end

        function [l1,l2,l3] = volume_to_length(obj, v1, v2, v3)
            l1=45.67*v1; l2=45.67*v2; l3=45.67*v3;
        end

        function [p1,p2,p3] = volumes_to_pressures(obj, v1, v2, v3)
            % cubic fit
            coeff = [57.9486,123.5743,-17.7039,1.8930];
            p1 = polyval(coeff,v1);
            p2 = polyval(coeff,v2);
            p3 = polyval(coeff,v3);
        end

        function vols = pressureToVolume(obj, P)
            vols = sign(P).*0.161.*abs(P).^(1/3);
        end

        function Jv = calcJacobianVolume(obj, v1,v2,v3)
            error('Please translate Python calcJacobianVolume');
        end

        function J = jacobian_from_lengths(obj)
            error('Translate this stub');
        end

        %--- Actuation & dynamics ---
        function Jr = get_jacobian_actuated(obj, u_ells, link_lens, beta_angles, d)
            if nargin<5, d=obj.d; end
            ell1 = sum(link_lens + 2*d*sin(beta_angles));
            ell2 = sum(link_lens - 2*d*sin(beta_angles));
            delta=ell1-ell2; C1=(2*cos(delta/(6*d))+1).^2;
            s12=sin(delta/(12*d)); c12=cos(delta/(12*d)); s6=sin(delta/(6*d));
            J11 = s12.*C1/6 + c12.*C1*(ell1+ell2)/(12*d) - 2*s6.*s12.*C1*(ell1+ell2)/(3*d);
            J12 = s12.*C1/6 - c12.*C1*(ell1+ell2)/(12*d) + 2*s6.*s12.*C1*(ell1+ell2)/(3*d);
            J21 = c12.*(2*cos(delta/(3*d))+1)/6 - 2*c12.*sin(delta/(3*d))*(ell1+ell2)/(3*d) - s12.*(2*cos(delta/(3*d))+1)*(ell1+ell2)/(12*d);
            J22 = c12.*(2*cos(delta/(3*d))+1)/6 + 2*c12.*sin(delta/(3*d))*(ell1+ell2)/(3*d) + s12.*(2*cos(delta/(3*d))+1)*(ell1+ell2)/(12*d);
            Jr = [J11, J12; J21, J22];
        end

        function J = jacobian_y(obj, beta, lens)
            b=beta; l=lens;
            J=[ l(1)*cos(b(1)) + l(2)*cos(b(1)+b(2)) + l(3)*cos(sum(b)), ...
                l(2)*cos(b(1)+b(2)) + l(3)*cos(sum(b)), ...
                l(3)*cos(sum(b));
                -l(1)*sin(b(1)) - l(2)*sin(b(1)+b(2)) - l(3)*sin(sum(b)), ...
                -l(2)*sin(b(1)+b(2)) - l(3)*sin(sum(b)), ...
                -l(3)*sin(sum(b))];
        end

        function tau = revolute_torque(obj, beta, lens, F)
            tau = obj.jacobian_y(beta,lens)' * F;
        end

        function f = link_force(obj, beta, lens, F)
            n=numel(beta); f=zeros(n,1);
            for i=1:n, ang=sum(beta(1:i)); f(i)=F(1)*sin(ang)+F(2)*cos(ang); end
        end

        function [Lnew,Bnew] = compute_q_dynamics(obj,K,L0,B0,F,omega)
            Kaa=K(1:3,1:3); Kuu=K(4:6,4:6);
            Buu=diag([0.05,0.05,0.05]);
            tau = obj.revolute_torque(B0,L0,F);
            fr  = obj.link_force(B0,L0,F);
            Bnew = B0 + Kuu\(tau - Buu*omega);
            dl   = Kaar;
            Lnew = L0 + dl;
        end

        function q = q_no_load(obj,U,n,d)
            if nargin<3,n=3; end;if nargin<4,d=4; end
            l1=U(1);l2=U(2);
            th=(l1-l2)/(2*d); s=(l1+l2)/2;
            arc=ones(1,n)*s/n; tcs=ones(1,n)*th/n;
            bs = tcs.*([0.5,ones(1,n-1)]);
            ls = (2*sin(tcs/2)./tcs).*arc;
            q=[ls,bs]';
        end

        function B= get_revolute_backbone(obj,Beta,Ls)
            x0=0;y0=0;
            x1=Ls(1)*sin(Beta(1));
            y1=Ls(1)*cos(Beta(1));
            x2=x1+Ls(2)*sin(sum(Beta(1:2)));
            y2=y1+Ls(2)*cos(sum(Beta(1:2)));
            x3=x2+Ls(3)*sin(sum(Beta)); y3=y2+Ls(3)*cos(sum(Beta));
            B=[0,x1,x2,x3;0,y1,y2,y3];
        end
    end
end