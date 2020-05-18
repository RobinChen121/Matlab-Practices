function ScenarioGeneration
% Management Science paper: generating scenario trees for multistage
% decision problems (2002)
 



iniDemand = [100, 200, 300, 500, 600];
iniPossibility = [0.2, 0.2, 0.2, 0.2, 0.2];

T = length(iniDemand);
% mean = 467.25;
% variance = 99422;
% worstDemand = mean - 2.5*sqrt(variance);
% Aeq =  [0, 0, 0, 0, 0, 1, 1, 1, 1, 1; 1, zeros(1, 2*T - 1); zeros(1, T), 1, zeros(1, T-1)];
% beq = [1; worstDemand; 0.01];
Aeq =  [0, 0, 0, 0, 0, 1, 1, 1, 1, 1,0, 0, 0, 0, 0,0, 0, 0, 0, 0];
beq = 1;
lb = ones(1, 4*T)*0.1;
% x = fmincon(@(x)ObjFun(x), [iniDemand, iniPossibility, iniDemand,iniDemand], [], [], Aeq, beq, lb, []);
[x, fval] = fmincon(@(x)ObjFun(x), [iniDemand, iniPossibility, iniDemand,iniDemand], [], [], Aeq, beq, lb, []);
d = x(1 : T);
p = x(T + 1 : 2*T);
d2 = x(2*T+1 : 3*T);
d3 = x(3*T+1 : 4*T);
fprintf('demand scenario:')
disp(d)
disp(d2)
disp(d3)
fprintf('possibility:')
disp(p)
disp(sum(p))
fprintf('fit gap = %.2f\n', fval);
fprintf('mean = %.2f\n', Mean(d, p));
fprintf('mean2 = %.2f\n', Mean(d2, p));
fprintf('mean3 = %.2f\n', Mean(d3, p));
fprintf('variance = %.2f\n', Variance(d, p));
fprintf('variance2 = %.2f\n', Variance(d2, p));
fprintf('variance3 = %.2f\n', Variance(d3, p));
fprintf('skew = %.2f\n', Skew(d, p));
fprintf('skew2 = %.2f\n', Skew(d2, p));
fprintf('skew3 = %.2f\n', Skew(d3, p));
fprintf('kurt = %.2f\n', Kurt(d, p));
fprintf('kurt2 = %.2f\n', Kurt(d2, p));
fprintf('kurt3 = %.2f\n', Kurt(d3, p));