function ScenarioGeneration
% Management Science paper: generating scenario trees for multistage
% decision problems (2002)
 



iniDemand = [400, 300, 100, 500, 600];
iniPossibility = [0.2, 0.2, 0.2, 0.2, 0.2];

T = length(iniDemand);
% mean = 467.25;
% variance = 99422;
% worstDemand = mean - 2.5*sqrt(variance);
% Aeq =  [0, 0, 0, 0, 0, 1, 1, 1, 1, 1; 1, zeros(1, 2*T - 1); zeros(1, T), 1, zeros(1, T-1)];
% beq = [1; worstDemand; 0.01];
Aeq =  [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
beq = 1;
lb = ones(1, 2*T)*0.01;
x = fmincon(@(x)ObjFun(x), [iniDemand, iniPossibility], [], [], Aeq, beq, lb, []);
[x, fval] = fmincon(@(x)ObjFun(x), [iniDemand, iniPossibility], [], [], Aeq, beq, ones(1, 2*T)*0.01, []);
d = x(1 : T);
p = x(T + 1 : 2*T);
fprintf('demand scenario:')
disp(d)
fprintf('possibility:')
disp(p)
fprintf('fit gap = %.2f\n', fval);
fprintf('mean = %.2f\n', Mean(d, p));
fprintf('variance = %.2f\n', Variance(d, p));
fprintf('skew = %.2f\n', Skew(d, p));
fprintf('kurt = %.2f\n', Kurt(d, p));