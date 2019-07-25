@@ -1,9 +1,17 @@
function obj = ObjFun(x)

mean = 467.25;
variance = 99.42;
skew = 1.06;
kurt = 4.35;
% mean = 467.25;
% variance = 99.42;
% skew = 1.06;
% kurt = 4.35;

mean = 33.82;
variance = 175.42;
skew = 0.25;
kurt = 2.78;


z

T = length(x)/2;
d = x(1:T);