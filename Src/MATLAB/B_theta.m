function [b,db_t, db_l] = B_theta(L,q2)

[f,d_f] = f_theta(q2);


if abs(q2) < 1e-3
    b = L/2;
    db_t = 0;  
    db_l = 1/2;

else
    b = (L/q2)*sqrt(f);
    db_t = (-L/(q2^2))*sqrt(f) + (L/(2*q2*sqrt(f)))*d_f;
    db_l = sqrt(f)/q2;
end