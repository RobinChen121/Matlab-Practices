function  RiskAversion
%RISKAVERSION 此处显示有关此函数的摘要
%   此处显示详细说明

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
con17 = (170*w(3,1) + 150 *w(3,2) + 36*w(3,3) + 10*w(3,4) - 170*1.4*y(3,1) - 150*1.4*y(3,2) - 150*x(1) - 230 *x(2) - 260*x(3) >= 48000);
constraints = [con1; con2; con3; con4; con5; con6; con7; con8; con9; con10; con11; con12; con13; con14; con15; con16; con17];

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
disp(108390 - value(obj))

end

