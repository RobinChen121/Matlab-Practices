function farmer
% /* total acres of land: 500
%  * minimum requirement: wheat >= 200 tons
%  * minimum requirement: corn >= 240 tons
%  * planning cost: wheat 150/ton, corn 230/ton, sugar 260/ton
%  * selling price: wheat 170/ton, corn 150/ton
%  * selling price: sugar 36/ton if yields <= 6000, 10/ton for excess  if yields > 6000
%  * yield rate: 2.5, 3, 20 respectively
%  * purchasing price: wheat 170*1.4/ton, corn 150*1.4/ton
%  * 
%  * decision variable: 
%  * x1: acres of land devoted to wheat
%  * x2: acres of land devoted to corn
%  * x3: acres of land devoted to sugar
%  * y1: tons of purchased wheat
%  * y2: tons of purchased corn
%  * w1: tons of selling wheat
%  * w2: tons of selling corn
%  * w3: tons of selling sugar at price 36
%  * w4: tons of selling sugar at price 10
%  * 
%  * build model:
%  * objective: max 170*w1 + 150 *w2 + 36*w3 + 10*w4 - 170*1.4*y1 - 150*1.4*y2 - 150*x1 - 230 *x2 - 260*x3
%  * 
%  * constraints:
%  * x1 + x2 + x3 <= 500
%  * 2.5*x1 + y1 - w1 >= 200
%  * 3*x2 + y2 -w2 >= 240
%  * w3 + w4 <= 20*x3
%  * w3 <= 6000
%  * x1, x2, x3, y1, y2, w1, w2, w3, w4 >= 0
%  * 
%  */

%% build yalmip model for normal yield
% decision variables
x = sdpvar(1,3);
y = sdpvar(1,2);
w = sdpvar(1,4);

% objective
obj = 170*w(1) + 150 *w(2) + 36*w(3) + 10*w(4) - 170*1.4*y(1) - 150*1.4*y(2) - 150*x(1) - 230 *x(2) - 260*x(3);

%constraints
con1 = (sum(x) <= 500);
con2 = (2.5*x(1) + y(1) - w(1) >= 200);
con3 = (3*x(2) + y(2) - w(2) >= 240);
con4 = (w(3) + w(4) <= 20*x(3));
con5 = (w(3) <= 6000);
con6 = (x >= 0 );
con7 = (y >= 0 );
con8 = (w >= 0);
constraints = [con1; con2; con3; con4; con5; con6; con7; con8];

% solve
diagnostics = optimize(constraints, -obj);
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



%% build yalmip model for yield scenarios
% 3 scenarios with equal possibility, high yield, normal yield, low yield
% decision variables
x = sdpvar(1,3);
y = sdpvar(3,2);
w = sdpvar(3,4);

% objective
obj = 1/3*(170*w(1,1) + 150 *w(1,2) + 36*w(1,3) + 10*w(1,4) - 170*1.4*y(1,1) - 150*1.4*y(1,2)) +...
         1/3*(170*w(2,1) + 150 *w(2,2) + 36*w(2,3) + 10*w(2,4) - 170*1.4*y(2,1) - 150*1.4*y(2,2)) +...
         1/3*(170*w(3,1) + 150 *w(3,2) + 36*w(3,3) + 10*w(3,4) - 170*1.4*y(3,1) - 150*1.4*y(3,2)) ...
         - 150*x(1) - 230 *x(2) - 260*x(3);

%constraints
con1 = (sum(x) <= 500);
con2 = (3*x(1) + y(1,1) - w(1,1) >= 200);
con3 = (3.6*x(2) + y(1,2) - w(1,2) >= 240);
con4 = (w(1,3) + w(1,4) <= 24*x(3));
con5 = (w(1,3) <= 6000);
con6 = (x >= 0 );
con7 = (y >= 0 );
con8 = (w >= 0);
con9 = (2.5*x(1) + y(2,1) - w(2,1) >= 200);
con10 = (2*x(1) + y(3,1) - w(3,1) >= 200);
con11 = (3*x(2) + y(2,2) - w(2,2) >= 240);
con12 = (2.4*x(2) + y(3,2) - w(3,2) >= 240);
con13 = (w(2,3) + w(2,4) <= 20*x(3));
con14 = (w(3,3) + w(3,4) <= 16*x(3));
con15 = (w(2,3) <= 6000);
con16 = (w(3,3) <= 6000);
constraints = [con1; con2; con3; con4; con5; con6; con7; con8; con9; con10; con11; con12; con13; con14; con15; con16];

% solve
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


%% price effect
% profit increase as price decrease
% decision variables
x = sdpvar(1,3);
y = sdpvar(1,2);
w = sdpvar(1,4);

% objective
obj = 170*1.1*w(1) + 150 *1.1*w(2) + 36*w(3) + 10*w(4) - 170*1.4*1.1*y(1) - 150*1.4*1.1*y(2) - 150*x(1) - 230 *x(2) - 260*x(3);

%constraints
con1 = (sum(x) <= 500);
con2 = (2.5*0.9*x(1) + y(1) - w(1) >= 200);
con3 = (3*0.9*x(2) + y(2) - w(2) >= 240);
con4 = (w(3) + w(4) <= 20*x(3));
con5 = (w(3) <= 6000);
con6 = (x >= 0 );
con7 = (y >= 0 );
con8 = (w >= 0);
constraints = [con1; con2; con3; con4; con5; con6; con7; con8];

% solve
diagnostics = optimize(constraints, -obj);
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