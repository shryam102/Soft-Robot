function xi = m_q_to_xi(phi, theta, L)
    eps = 1e-6;
    if abs(theta) < eps
        b = L / 2;
        eta = 0;
    else
        s = sin(theta / 2);
        b = (L / theta) * sqrt(1 + 4 * (s / theta)* ((s / theta) - cos(theta / 2)));
        eta = acos((L / theta) * s / b);
    end

    xi = [
        phi;
        theta / 2 - eta;
        b;
        eta;
        -phi;
        phi;
        eta;
        b;
        theta / 2 - eta;
        -phi
    ];