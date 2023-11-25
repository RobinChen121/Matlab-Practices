function RestrictLineConvex

figure(1);
fplot(@(x)x.^2, [-2 2]);

hold on;
x0 = 0;
v = 2;
fplot(@(t) (x0 + t*v).^2, [-2 2]);
legend('f(x)', 'f(t)=f(x+tv)', 'Location', 'north')

figure(2);
f = @(x, y) x.^2 + y.^2;
fmesh(f)
xlabel('x');
ylabel('y');
zlabel('z');

hold on;
x0 = [0  0];
v = [1 1];
fmesh(@(t) t);


end



