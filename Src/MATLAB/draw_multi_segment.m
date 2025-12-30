function draw_multi_segment(q_list, curves,color)
    L = curves.L;
    n = numel(q_list)/2;
    T = eye(4);
    all_points = zeros(3, 1);
    thetas = q_list(2:2:2*n);
    phis = q_list(1:2:2*n);
    ls = linspace(L,L,length(thetas));
%     l_cc = L*n;
%     theta_cc = q_cc(2);
%     phi_cc = q_cc(1);

    for i = 1:n
        phi = q_list(2*i-1);
        theta = q_list(2*i);
        xi = m_q_to_xi(phi, theta, L);

        dh_params = [
            xi(1),    0,     0, -pi/2;
            xi(2),    0,     0,  pi/2;
            0,        xi(3), 0, -pi/2;
            xi(4),    0,     0,  pi/2;
            xi(5),    0,     0,     0;
            xi(6),    0,     0, -pi/2;
            xi(7),    0,     0,  pi/2;
            0,        xi(8), 0, -pi/2;
            xi(9),    0,     0,  pi/2;
            xi(10),   0,     0,     0;
        ];

        for j = 1:10
            th = dh_params(j,1);
            d  = dh_params(j,2);

            a  = dh_params(j,3);
            al = dh_params(j,4);

            Tj = DH_transform(th, d, a, al);
            T = T * Tj;
            all_points(:, end+1) = T(1:3,4);
        end
    end

    plot3(all_points(1,:), all_points(2,:), all_points(3,:), '--o','LineWidth', 1.5,'Color', [0 0 0 0.3] );
    hold on;
    FW_Curve(thetas,phis,ls, color);
%     FW_Curve(theta_cc,phi_cc,l_cc,'b');
%     legend({'Augemented Model','PCC','Single Segment'});
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    xlim([-L*n,L*n]);
    ylim([-L*n,L*n]);
    zlim([0,L*n]);
    grid on;
    drawnow; 
    hold off
