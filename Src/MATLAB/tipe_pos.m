function r = tipe_pos(q)
n = numel(q)/3;
r = zeros(3,1);
for i = 1:n
    idx  = 3*(i-1) + (1:3);
    phi = q(idx(1));
    theta = q(idx(2));
    s = q(idx(3));

    if abs(theta) < 1e-3
        ri = [0;0;s];
    else
        ri = [s*cos(phi)*(1 - cos(theta))/theta;
              s*sin(phi)*(1 - cos(theta))/theta;
              s*sin(theta)/theta];
    end
    r = r + ri;
end
end