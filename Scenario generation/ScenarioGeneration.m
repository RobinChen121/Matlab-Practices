function ScenarioGeneration
% this a matlab test code for getting the demand realizaions in the scenario tree. This 
% method is adopted from the paper in Hoyland and Wallace (2001),
% Management Science
%   

% mean = 467.25
% variance = 99.42
% skewness = 1.06
% kurosis = 4.35

iniDemand = [400, 300, 100, 500, 600];
iniPossibility = [0.2, 0.2, 0.2, 0.2, 0.2];

T = length(iniDemand);
Aeq =  [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
beq = 1;
x = fmincon(@(x)ObjFun(x), [iniDemand, iniPossibility], [], [], Aeq, beq, zeros(1, 2*T),[]);
d = x(1 : T);
p = x(T + 1 : 2*T);
disp(d);
disp(p);
end

