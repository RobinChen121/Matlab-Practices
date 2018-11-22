function DrawHyperPlane
fmesh(@(x,y)1/4*(9-2*x-3*y));
hold on;

plot3([-1,3], [-2, 4], [-3, 5]);

title('Hyper plane: 2x+3y+4z=9')
str = '$$\frac{x-1}{2}=\frac{y-1}{3}=\frac{z-1}{4}$$';
text(-0.7, 9.8, str, 'Interpreter','latex');
annotation('arrow','X',[0.4,0.5],'Y',[0.6,0.6]

end