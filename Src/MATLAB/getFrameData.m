function [skelPts, softPts] = getFrameData(q_list, curves)
    L = curves.L;
    n = numel(q_list)/2;
    T = eye(4);
    skelPts = zeros(3,1);

    thetas = q_list(2:2:end);
    phis   = q_list(1:2:end);
    ls     = repmat(L, 1, n);

    % build skeleton points
    for i = 1:n
        xi = m_q_to_xi(phis(i), thetas(i), L);
        dh = [ xi(1),0,0,-pi/2; xi(2),0,0, pi/2;
               0, xi(3),0,-pi/2; xi(4),0,0, pi/2;
               xi(5),0,0,    0; xi(6),0,0,-pi/2;
               xi(7),0,0, pi/2; 0, xi(8),0,-pi/2;
               xi(9),0,0, pi/2; xi(10),0,0,   0 ];
        for j = 1:size(dh,1)
            T = T * DH_transform(dh(j,1), dh(j,2), dh(j,3), dh(j,4));
            skelPts(:,end+1) = T(1:3,4);
        end
    end

    % build soft‐curve points
    % (assuming you can adapt FW_Curve to return data instead of plotting)
    softPts = FW_Curve(thetas, phis, ls, 'collect');  
    % — or, if FW_Curve doesn’t support that, copy its internal logic
    %   here so you return a 3×M array of curve points.
end
