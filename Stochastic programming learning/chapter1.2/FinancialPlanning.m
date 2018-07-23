function FinancialPlanning

% initial capital b=55000
% goal final capital > G = 80000
% 3 periods
% 2 scenarios in each period, 
% scenario 1: profit rate 1.25 in stocks and 1.14 in bonds
% scenario 2: profit rate 1.06 in stocks and 1.12 in bonds
% there are 8 scenarios in total in 3 periods

% decision variables
% excess final capital y( s1, s2, s3)
% shortage final capital w(s1, s2, s3)
% amount of capital to item i in period t for scenario St-1: x(i, t, St-1)

% constraints
% intial capital sum; sum of initial scenarios equal to latter sum
% scenarios; final capital sum

%% build yalmip model
% variables
y = sdpvar(2, 2, 2);
w = sdpvar(2, 2, 2);
x1 = sdpvar(1, 2);
x2 = sdpvar(2, 2);  % xt(i, j) t½×¶Î item i scenario Sj
x3 = sdpvar(2, 2, 2);

% objective
obj = sum(sum(sum(0.125*y-4*w)));

% constraints
con1 = (sum(x1) == 55000);
con2 = (1.25*x1(1) + 1.14*x1(2) == x2(1, 1) +x2(2, 1));
con3 = (1.06*x1(1) + 1.12*x1(2) == x2(1,2) + x2(2, 2));
con4 = (1.25*x2(1, 1) + 1.14*x2(2, 1) == x3(1, 1, 1) +x3(2, 1, 1));
con5 = (1.06*x2(1, 1) + 1.12*x2(2, 1) == x3(1, 1, 2) + x3(2, 1, 2));
con6 = (1.25*x2(1, 2) + 1.14*x2(2, 2) == x3(1, 2, 1) +x3(2, 2, 1));
con7 = (1.06*x2(1, 2) + 1.12*x2(2, 2) == x3(1, 2, 2) + x3(2, 2, 2));
con8 = (1.25*x3(1, 1, 1) + 1.14*x3(2, 1, 1)- y(1, 1, 1) + w(1, 1, 1) == 80000);
con9 = (1.06*x3(1, 1, 1) + 1.12*x3(2, 1, 1)- y(1, 1, 2) + w(1, 1, 2) == 80000);
con10 = (1.25*x3(1, 1, 2) + 1.14*x3(2, 1, 2)- y(1, 2, 1) + w(1, 2, 1) == 80000);
con11 = (1.06*x3(1, 1, 2) + 1.12*x3(2, 1, 2)- y(1, 2, 2) + w(1, 2, 2) == 80000);
con12 = (1.25*x3(1, 2, 1) + 1.14*x3(2, 2, 1)- y(2, 1, 1) + w(2, 1, 1) == 80000);
con13 = (1.06*x3(1, 2, 1) + 1.12*x3(2, 2, 1)- y(2, 1, 2) + w(2, 1, 2) == 80000);
con14 = (1.25*x3(1, 2, 2) + 1.14*x3(2, 2, 2)- y(2, 2, 1) + w(2, 2, 1) == 80000);
con15 = (1.06*x3(1, 2, 2) + 1.12*x3(2, 2, 2)- y(2, 2, 2) + w(2, 2, 2) == 80000);
con16 = (y >= 0);
con17 = (w >= 0);
con18 = (x1 >= 0);
con19 = (x2 >= 0);
con20 = (x3 >= 0);

constraints = [con1; con2; con3; con4; con5; con6; con7; con8; con9; con10; con11; con12; con13; con14; con15; con16; con17; con18; con19; con20];

% solve
options = sdpsettings('solver', 'cplex');
diag = optimize(constraints, -obj, options);
if diag.problem == 0
    disp('Solver thinks it is feasible')
    value(obj)
    x1 = value(x)
    x2 = value(x2)
    x3 = value(x3)
    y = value(y)
    w = value(w)
elseif diag.problem == 1
    disp('Solver thinks it is infeasible')
else
    disp('Something else happened')
end

end