function  DrawContour

load GA

Y = 15:50;
X = 0:50;
contour(X, Y, GA, 10)
xlabel('y', 'FontWeight', 'bold');
ylabel('R', 'FontWeight', 'bold');
title('GA')
figure
mesh(X, Y, GA)
end

