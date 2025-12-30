function R = angle_axis(phi, theta)
    
    R = Rz(phi)*Ry(theta)*Rz(-phi);


end