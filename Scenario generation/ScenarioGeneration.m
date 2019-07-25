function ScenarioGeneration
% this a matlab test code for getting the demand realizaions in the scenario tree. This 
% method is adopted from the paper in Hoyland and Wallace (2001),
% Management Science
%   

% mean = 467.25;
% variance = 99.42;
mean = 33.82;
variance = 175.42;

<<<<<<< HEAD
iniDemand = [50, 300, 100, 80, 600];
=======
iniDemand = ones(1, 5) * mean/5;
>>>>>>> 9b61c7797d5db529eac904f64871378c22468c06
iniPossibility = [0.2, 0.2, 0.2, 0.2, 0.2];

T = length(iniDemand);
worstDemand = mean - 2.5*sqrt(variance);
% Aeq =  [0, 0, 0, 0, 0, 1, 1, 1, 1, 1; 1, zeros(1, 2*T - 1); zeros(1, T), 1, zeros(1, T-1)];
% beq = [1; worstDemand; 0.01];
Aeq =  [0, 0, 0, 0, 0, 1, 1, 1, 1, 1];
beq = 1;
<<<<<<< HEAD
ub = [10000*ones(1,T), 0.5*ones(1, T)];
[x, fval] = fmincon(@(x)ObjFun(x), [iniDemand, iniPossibility], [], [], Aeq, beq, ones(1, 2*T)*0.01, ub);
=======
[x, fval] = fmincon(@(x)ObjFun(x), [iniDemand, iniPossibility], [], [], Aeq, beq, ones(1, 2*T)*0.01, []);
>>>>>>> 9b61c7797d5db529eac904f64871378c22468c06
d = x(1 : T);
p = x(T + 1 : 2*T);
d
p
<<<<<<< HEAD
fprintf('fit value = %.2f\n', fval);
=======
fprintf('fit gap = %.2f\n', fval);
>>>>>>> 9b61c7797d5db529eac904f64871378c22468c06
fprintf('mean = %.2f\n', Mean(d, p));
fprintf('variance = %.2f\n', Variance(d, p));
fprintf('skew = %.2f\n', Skew(d, p));
fprintf('kurt = %.2f\n', Kurt(d, p));
end

