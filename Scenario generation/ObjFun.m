function obj = ObjFun(x)


% mean = 33.82;
% variance = 175.42;
% skew = 0.25;
% kurt = 2.78;

mean = [82.14, 397.33, 181.54];
variance = [5173.98, 79206.22, 22520.71];
skew = [2.06, 0.41, 1.67];
kurt = [7.39, 1.78, 5.37];

T = length(x)/4;
d = x(1:T);
p = x(T+1 : 2*T);
d2 = x(2*T+1 : 3*T);
d3 = x(3*T+1 : 4*T);

meanE = Mean(d, p);
varianceE = Variance(d, p);
skewE = Skew(d, p);
kurtE = Kurt(d, p);
meanE2 = Mean(d2, p);
varianceE2 = Variance(d2, p);
skewE2 = Skew(d2, p);
kurtE2 = Kurt(d2, p);
meanE3 = Mean(d3, p);
varianceE3 = Variance(d3, p);
skewE3 = Skew(d3, p);
kurtE3 = Kurt(d3, p);

obj = 0.1*(meanE - mean(1))^2 + 0.0001*(varianceE - variance(1))^2 + 4*(skewE - skew(1))^2 + 4*(kurtE - kurt(1))^2 + ...
    0.1*(meanE2 - mean(2))^2 + 0.0001*(varianceE2 - variance(2))^2 + 4*(skewE2 - skew(2))^2 + 4*(kurtE2 - kurt(2))^2 + ...
    0.1*(meanE3 - mean(3))^2 + 0.0001*(varianceE3 - variance(3))^2 + 4*(skewE3 - skew(3))^2 + 4*(kurtE3 - kurt(3))^2;


end
