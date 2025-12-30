% A Function which outputs the points in 3D space to reconstruct the SBA 
% based on the configuration variables or arc parameters (phi, theta, S)

function [softPts, skelPts, r] = compute_Soft_curve(q)
n = numel(q)/3;                  % n is the number of segments our system is currently divided into 
samples = 100;                   % number of points/ coordinates to be plotted for each segment
T = eye(4);                 
softPts = zeros(3,samples*n);    % Matrix to store points/coordinates belonging to soft robot (basically the curve)
skelPts = zeros(3, n+1);         % Matrix to store points/coordinates belonging to rigid robot 
phis = q(1:3:end);
thetas = q(2:3:end);
L = q(3:3:end);

for i =1:n
    Xi = position_vector(thetas(i), samples, phis(i), L(i), T);
    idx = (i-1)*samples + (1:samples);
    softPts(:,idx) = Xi;
    T = T*Transform(thetas(i), phis(i), L(i));
    skelPts(:, i+1) = softPts(:,samples*i);

end
r = T;              % Its the tip of the robot
end