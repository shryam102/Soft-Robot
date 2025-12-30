function X = position_vector(theta,n,phi,l,T)

t = linspace(0,theta,n);
X = zeros(3,n);
L = linspace(0,l,n);

for i = 1:n
    Mat = Transform(t(i), phi, L(i));

    Tf = T*Mat;

    X(:, i) = Tf(1:3, 4);
end

