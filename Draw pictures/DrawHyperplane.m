function DrawHyperplane

N = 50;
y = linspace(0,20,N);
x = linspace(0,20,N);
z = linspace(0,50,N);
[xx,yy] = meshgrid(x,y);
zz = 4 - xx - yy;
figure
surf(xx, yy, zz)
grid on

end