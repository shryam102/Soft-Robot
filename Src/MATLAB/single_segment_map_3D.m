function [xi, Jmi] = single_segment_map_3D(L, qi)

    q1 = qi(1);
    q2 = qi(2);
    [N,d_N,~] = N_theta(L,q2);
    [B,d_B,~] = B_theta(L,q2);

    xi = [      q1;
           (q2/2) - N;
                 B;
                 N;
               -q1;
                q1;
                 N;
                 B;
            (q2/2) - N;
                 -q1];

    Jmi = [1        0;
           0   ((1/2) - d_N);
           0       d_B;
           0       d_N;
          -1        0;
           1        0;
           0       d_N;
           0       d_B;
           0   ((1/2) - d_N);
          -1        0];

end