function [xi, Jm] = build_XI_J(q)
n = numel(q)/3;
xi = zeros(5*n,1);
Jm = zeros(5*n,3*n);
for i = 1:n
    idx = (i-1)*5 + (1:5);
    [xi_block, Jmi] = single_map(q(3*i-2:3*i));
    xi(idx) = xi_block;
    Jm(idx, 3*i-2:3*i) = Jmi;
end

end