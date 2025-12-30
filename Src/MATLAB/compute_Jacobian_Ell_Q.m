function J = compute_Jacobian_Ell_Q(q)

n = numel(q)/3;

J = zeros(3,3*n);
d = 4;

for i = 1:3
    for j = 1:3:3*n
        J(i,j) = compute_dEll_dphi(q,d,i,j);
    end
    for j = 2:3:3*n
        J(i,j) = compute_dEll_dtheta(q,d,i,j);
    end
    for j = 3:3:3*n
        J(i,j) = 1;
    end
end

    function x =compute_dEll_dphi(q,d,i,j)
        
        if i ==1
            x = -d*q(j+1)*cos(q(j));
        elseif i==2
            x =  d*q(j+1)*sin(pi/6 - q(j));
        else
            x =  d*q(j+1)*sin(pi/6 + q(j));
        end
    end

    function x = compute_dEll_dtheta(q,d,i,j)
       
        if i ==1
            x = -d*sin(q(j-1));
        elseif i ==2
            x =  d*cos(pi/6 - q(j-1));
        else
            x = -d*cos(pi/6 + q(j-1));
        end

    end

end