function CapacityExpansionStochastic
%% stochastic

% 使用情景方法处理了随机变量

% parameters
sd = [3, 5, 7];   % random demands
pd = [0.3, 0.4, 0.3];   % demand possiblities
tao = [10, 6, 1]; % load durations

% decision variables
x = sdpvar(1, 4);
y = sdpvar(4, 3, 3);

% objective function
obj1 = 10*x(1) + 7*x(2) + 16*x(3) + 6*x(4);
obj2 = 0;
for j = 1:3
    for s = 1:3
        obj2 = obj2 + pd(s)*tao(j)*(4*y(1,j,s) + 4.5*y(2,j,s) +3.2*y(3,j,s) + 5.5*y(4,j,s));
    end
end
obj = obj1 + obj2;


% constraints
con1 = (10*x(1) + 7*x(2) + 16*x(3) + 6*x(4) <= 120);
con21 =(sum(y(1,:,:)) <= x(1)) ;
con22 =(sum(y(2,:,:)) <= x(2)) ;
con23 =(sum(y(3,:,:)) <= x(3)) ;
con24 =(sum(y(4,:,:)) <= x(4)) ;
con31 = (sum(y(:,1,1)) == sd(1));
con32 = (sum(y(:,1,2)) == sd(2));
con33 = (sum(y(:,1,3)) == sd(3));
con41 = (sum(y(:, 2, :)) == 3);
con42 = (sum(y(:, 3, :)) == 2);
con5 = (x >= 0);
con6 = (y >= 0);
con7 = (sum(x) >=12);

constraints = [con1; con21; con22; con23; con24; con31; con32; con33; con41; con42; con5; con6; con7];

% solve
diagnostics = optimize(constraints, obj);
if diagnostics.problem == 0
    disp('Solver thinks it is feasible')
    value(obj)
    x = value(x)
    y = value(y)
elseif diagnostics.problem == 1
    disp('Solver thinks it is infeasible')
else
    disp('Something else happened')
end
end