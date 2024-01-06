function RestrictLineConvex

close all;
figure(1);
fplot(@(x)x.^2, [-2 2]);

hold on;
x0 = 0;
v = 2;
fplot(@(t) (x0 + t*v).^2, [-2 2]);
v = 0.5;
fplot(@(t) (x0 + t*v).^2, [-2 2]);
x0 = 1;
v = 2;
fplot(@(t) (x0 + t*v).^2, [-2 2]);
legend('f(x)=x^2', 'f(t)=f(x+tv)=4t^2', 'f(t)=f(x+tv)=0.25t^2', 'f(t)=f(x+tv)=(1+2t)^2', 'Location', 'north')

figure(2);
f = @(x, y) x.^2 + y.^2;
fmesh(f)
xlabel('x');
ylabel('y');
zlabel('z');
title('$f(x) = x_1^2+x_2^2$', 'interpreter', 'latex')

hold on;
x0 = [0  0];
v = [1 1];
patch([-5 -5 5 5], [-5 -5 5 5], [-10 50 50 10], 'r')

figure(3)
fplot(@(t) (2*t).^2, [-2 2]);
title('$f(t)=f(x+tv)=2t^2$', 'interpreter', 'latex')
end



