function DrawSaddlePoint

% draw saddle points

% f(x)=x^3
x=-2:0.1:2;
y=x.^3;
plot(x,y);
axis equal;
hold on;
scatter(0,0,'filled','MarkerFaceColor','r');
title('f(x)=x^3');

% f(x)=x^2-y^2
figure (2);
[x,y]=meshgrid(-2:0.03:2);
z=x.^2-y.^2;
mesh(x,y,z);
hold on;
scatter3(-2,-2,0,'filled','MarkerFaceColor','r');
scatter3(-1,-1,0,'filled','MarkerFaceColor','r');
scatter3(0,0,0,'filled','MarkerFaceColor','r');
scatter3(1,1,0,'filled','MarkerFaceColor','r');
scatter3(2,2,0,'filled','MarkerFaceColor','r');
text(0,0,0,'\leftarrow');
xlabel('x');
ylabel('y');
zlabel('z');
title('f(x,y)=x^2-y^2');
hold off;
end