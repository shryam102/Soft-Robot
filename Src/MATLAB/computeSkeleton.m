function skelPts = computeSkeleton(q, curves)
    L = curves.L;
    n = numel(curves);
    T = eye(4);
    pointsPerSegment = 10;               
    totalPts = n * pointsPerSegment;
    skelPts = zeros(3, totalPts);

    k = 1;

    for i = 1:n
        phi   = q(2*i-1);
        theta = q(2*i);
        xi    = m_q_to_xi(phi, theta, L);

        dh = [ xi(1),   0,     0, -pi/2; 
               xi(2),   0,     0,  pi/2; 
               0,       xi(3), 0, -pi/2; 
               xi(4),   0,     0,  pi/2;
               xi(5),   0,     0,     0; 
               xi(6),   0,     0, -pi/2;
               xi(7),   0,     0,  pi/2; 
               0,      xi(8),  0, -pi/2;
               xi(9),   0,     0,  pi/2; 
               xi(10),  0,     0,    0 ];

        for j = 1:size(dh,1)
            T = T * DH_transform(dh(j,1), dh(j,2), dh(j,3), dh(j,4));
            skelPts(:,k) = T(1:3,4);
            k= k + 1;
        end
    end
end
