function ScenarioGeneration
% this a matlab test code for getting the demand realizaions in the scenario tree. This 
% method is adopted from the paper in Hoyland and Wallace (2001),
% Management Science
%   

% mean = 467.25
% variance = 99.42
% skewness = 1.06
% kurosis = 4.35

iniDemand = [50, 300, 100, 80, 600];
iniPossibility = [0.2, 0.2, 0.2, 0.2, 0.2];

T = length(iniDemand);
Aeq =  [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
beq = 1;
ub = [10000*ones(1,T), 0.5*ones(1, T)];
[x, fval] = fmincon(@(x)ObjFun(x), [iniDemand, iniPossibility], [], [], Aeq, beq, ones(1, 2*T)*0.01, ub);
d = x(1 : T);
p = x(T + 1 : 2*T);
d
p
fprintf('fit value = %.2f\n', fval);
fprintf('mean = %.2f\n', Mean(d, p));
fprintf('variance = %.2f\n', Variance(d, p));
fprintf('skew = %.2f\n', Skew(d, p));
fprintf('kurt = %.2f\n', Kurt(d, p));
end

