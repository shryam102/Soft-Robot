function links = dh_params()
dh = [ 
     pi/4,     0.2,   0,   -pi/2,   1.0;  
     0,        0.3,   0,    pi/2,   2.5;  
    -pi/6,     0.0,   0,     0,     1.2; 
     ];

n = size(dh,1);
links(n) = struct();  
for i = 1:n
    links(i).theta = dh(i,1);
    links(i).d     = dh(i,2);
    links(i).a     = dh(i,3);
    links(i).alpha = dh(i,4);
    links(i).mu    = dh(i,5);
end
end
