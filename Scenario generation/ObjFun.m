function obj = ObjFun(x)


% mean = 33.82;
% variance = 175.42;
% skew = 0.25;
% kurt = 2.78;

mean = 467.25;
variance = 99422;
skew = 1.06;
kurt = 4.35;

T = length(x)/2;
d = x(1:T);
p = x(T+1 : 2*T);
meanE = Mean(d, p);
varianceE = Variance(d, p);
skewE = Skew(d, p);
kurtE = Kurt(d, p);

obj = 0.1*(meanE - mean)^2 + 0.1*(varianceE - variance)^2 + 0.4*(skewE - skew)^2 ...
        + 0.4*(kurtE - kurt)^2;


end
