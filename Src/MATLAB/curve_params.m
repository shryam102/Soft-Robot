function curves = curve_params()
n         = 3;      
segLength = 5;  
segMass   = 0.20; 

curves(n) = struct('L',[],'mu',[]);
for i = 1:n
    curves(i).L  = segLength;
    curves(i).mu = segMass;
end
end
