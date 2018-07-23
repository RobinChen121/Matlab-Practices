function  DataFluc
%DATAFLUC 此处显示有关此函数的摘要
%   purchase price of wheat being 300 or 220 with equal probability

%  * decision variable: 
%  * x1: acres of land devoted to wheat
%  * x2: acres of land devoted to corn
%  * x3: acres of land devoted to sugar
%  * y1: tons of purchased wheat
%  * y2: tons of purchased corn
%  * w1: tons of selling wheat
%  * w2: tons of selling corn
%  * w3: tons of selling sugar at price 41
%  * w4: tons of selling sugar at price 11

%% build yalmip model for yield scenarios
% decision variables
x = sdpvar(1,3);
y = sdpvar(2,2);
w = sdpvar(2,4);

% objective
obj = 1/2*(300*w(1,1) + 170 *w(1,2) + 41*w(1,3) + 11*w(1,4) - 300*1.4*y(1,1) - 170*1.4*y(1,2)) +...
         1/2*(220*w(2,1) + 170 *w(2,2) + 41*w(2,3) + 11*w(2,4) - 220*1.4*y(2,1) - 170*1.4*y(2,2)) +...
         - 180*x(1) - 280 *x(2) - 310*x(3);

%constraints
con1 = (sum(x) <= 500);
con2 = (2.5*x(1) + y(1,1) - w(1,1) >= 200);
con3 = (2.5*x(1) + y(2,1) - w(2,1) >= 200);
con4 = (3*x(2) + y(1,2) - w(1,2) >= 240);
con5 = (3*x(2) + y(2,2) - w(2,2) >= 240);
con6 = (w(2,3) + w(2,4) <= 20*x(3));
con7 = (w(1,3) + w(1,4) <= 20*x(3));
con8 = (w(1,3) <= 6000);
con9 = (w(2,3) <= 6000);
con10 = (x >= 0 );
con11 = (y >= 0 );
con12 = (w >= 0);
constraints = [con1; con2; con3; con4; con5; con6; con7; con8; con9; con10; con11; con12];

% solve
%options = sdpsettings('solver','CPLEX');
diagnostics = optimize(constraints, -obj);
value(obj)
if diagnostics.problem == 0
    disp('Solver thinks it is feasible')
    value(obj)
    x = value(x)
    y = value(y)
    w = value(w)
elseif diagnostics.problem == 1
    disp('Solver thinks it is infeasible')
else
    disp('Something else happened')
end

y = sdpvar(2,2);
w = sdpvar(2,4);

% objective
obj = 1/2*(300*w(1,1) + 170 *w(1,2) + 41*w(1,3) + 11*w(1,4) - 300*1.4*y(1,1) - 170*1.4*y(1,2)) +...
         1/2*(220*w(2,1) + 170 *w(2,2) + 41*w(2,3) + 11*w(2,4) - 220*1.4*y(2,1) - 170*1.4*y(2,2)) +...
         - 180*x(1) - 280 *x(2) - 310*x(3);

%constraints
con1 = (sum(x) <= 500);
con2 = (2.5*x(1) + y(1,1) - w(1,1) >= 200);
con3 = (2.5*x(1) + y(2,1) - w(2,1) >= 200);
con4 = (3*x(2) + y(1,2) - w(1,2) >= 240);
con5 = (3*x(2) + y(2,2) - w(2,2) >= 240);
con6 = (w(2,3) + w(2,4) <= 20*x(3));
con7 = (w(1,3) + w(1,4) <= 20*x(3));
con8 = (w(1,3) <= 6000);
con9 = (w(2,3) <= 6000);
con10 = (x >= 0 );
con11 = (y >= 0 );
con12 = (w >= 0);
constraints = [con1; con2; con3; con4; con5; con6; con7; con8; con9; con10; con11; con12];

% solve
%options = sdpsettings('solver','CPLEX');
diagnostics = optimize(constraints, -obj);
value(obj)
if diagnostics.problem == 0
    disp('Solver thinks it is feasible')
    value(obj)
    x = value(x)
    y = value(y)
    w = value(w)
elseif diagnostics.problem == 1
    disp('Solver thinks it is infeasible')
else
    disp('Something else happened')
end
end

