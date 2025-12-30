% Function which outputs the transformation matrix based on the theta, phi
% and arc length l (or S)

% The homogenous trasnformation for a single segment having theta, phi, l
% can be evaluated considering three coordinate transformation
% First, rotation about z axis by phi degrees
% Second, roatation of theta about y axis and then moving to the end point
% of the robot
% We could have left it here, but its useful to orient the tip frame such
% that it alligns with the base frame when it 'slids' along the curve
% without rotation about the z axis. To achieve this we simply re-rotate
% by -phi about the z axis.
% Thus our transformation matrix is achieved


function x = Transform(theta,phi,l)

if theta == 0
    x = [Rz(phi), [0;0;0]; 0 0 0 1]*[Ry(theta), l*[0; 0; 1]; 0 0 0 1]*[Rz(-phi), [0;0;0]; 0 0 0 1];
else
    r = l/theta;

    x = [Rz(phi), [0;0;0]; 0 0 0 1]*[Ry(theta), r*[(1 - cos(theta)); 0; sin(theta)]; 0 0 0 1]*[Rz(-phi), [0;0;0]; 0 0 0 1];

end
