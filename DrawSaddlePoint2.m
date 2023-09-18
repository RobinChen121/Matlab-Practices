function DrawSaddlePoint2

h = 1;
p = 2;
K = 20;
lambda = 10;
x = linspace(0, 1, 20);
q = linspace(1, 200, 20);

[X,Q] = meshgrid(x,q);
Z = h.*Q.*(1 - X.*X)/2 + p*Q.*X.*X/2 + K*lambda./Q;
Z = Q.*X.^2;
mesh(X,Q,Z);

xlabel('x');
ylabel('Q');
zlabel('g');
end

