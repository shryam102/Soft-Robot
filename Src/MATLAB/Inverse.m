function i_k = Inverse(X)

phi = atan2(X(2), X(1));

x_dash = [X(1); X(2)];

k = 2*norm(x_dash)/norm(X)^2;


if X(3) > 0
    
    theta = acos(1 - norm(x_dash)*k);

elseif X(3) <= 0

    theta = 2*pi - acos(1 - norm(x_dash)*k);
end

s = theta*k ;
    
if X(1) == 0 && X(2) == 0
    phi = 0;
    k = 0;
    s = X(3);
end

i_k = [k; theta; phi; s];






