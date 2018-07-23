function BinaryFirstStage

%% problem description
% similar to the farmer's problem, but 4 fixed size fields
% and each field raises only one crop

%  minimum requirement: wheat >= 200 tons
%  minimum requirement: corn >= 240 tons
%  planning cost: wheat 150/ton, corn 230/ton, sugar 260/ton
%  selling price: wheat 170/ton, corn 150/ton
%  selling price: sugar 36/ton if yields <= 6000, 10/ton for excess  if yields > 6000
%  yield rate: 2.5, 3, 20 respectively at normal yield scenario
%  purchasing price: wheat 170*1.4/ton, corn 150*1.4/ton  

%  * decision variable: 
%  * x11: grow wheat in land 1
%  * x12: grow corn in land 1
%  * x13: grow sugar in land 1
%  * x21: grow wheat in land 2
%  * x22: grow corn in land 2
%  * x23: grow sugar in land 2
%  * x31: grow wheat in land 3
%  * x32: grow corn in land 3
%  * x33: grow sugar in land 3
%  * x41: grow wheat in land 4
%  * x42: grow corn in land 4
%  * x43: grow sugar in land 4

%  * y1s: tons of purchased wheat at scenario s
%  * y2s: tons of purchased corn at scenario s
%  * w1s: tons of selling wheat at scenario s
%  * w2s: tons of selling corn at scenario s
%  * w3s: tons of selling sugar at price 36 at scenario s
%  * w4s: tons of selling sugar at price 10 at scenario s

%% build yalmip model for 3 yield scenarios
% 3 scenarios with equal possibility, high yield, normal yield, low yield
% decision variables
x = binvar(4,3);
y = sdpvar(3,2);
w = sdpvar(3,4);

% objective
obj = 1/3*(170*w(1,1) + 150 *w(1,2) + 36*w(1,3) + 10*w(1,4) - 170*1.4*y(1,1) - 150*1.4*y(1,2)) +...
         1/3*(170*w(2,1) + 150 *w(2,2) + 36*w(2,3) + 10*w(2,4) - 170*1.4*y(2,1) - 150*1.4*y(2,2)) +...
         1/3*(170*w(3,1) + 150 *w(3,2) + 36*w(3,3) + 10*w(3,4) - 170*1.4*y(3,1) - 150*1.4*y(3,2)) ...
         - 150*(185*x(1,1)+145*x(2,1)+105*x(3,1)+65*x(4,1)) -...
         230 *(185*x(1,2)+145*x(2,2)+105*x(3,2)+65*x(4,2)) -...
         260*(185*x(1,3)+145*x(2,3)+105*x(3,3)+65*x(4,3));
 
% constraints
con1 = (sum(x(1,:)) == 1);
con2 = (sum(x(2,:)) == 1);
con3 = (sum(x(3,:)) == 1);
con4 = (sum(x(4,:)) == 1);
con5 = (3*(185*x(1,1)+145*x(2,1)+105*x(3,1)+65*x(4,1)) + y(1,1) - w(1,1) >= 200);
con6 = (3.6*(185*x(1,2)+145*x(2,2)+105*x(3,2)+65*x(4,2)) + y(1,2) - w(1,2) >= 240);
con7 = (w(1,3) + w(1,4) <= 24*(185*x(1,3)+145*x(2,3)+105*x(3,3)+65*x(4,3)));
con8 = (w(1,3) <= 6000);
con9 = (y >= 0 );
con10 = (w >= 0);
con11 = (2.5*(185*x(1,1)+145*x(2,1)+105*x(3,1)+65*x(4,1)) + y(2,1) - w(2,1) >= 200);
con12 = (2*(185*x(1,1)+145*x(2,1)+105*x(3,1)+65*x(4,1)) + y(3,1) - w(3,1) >= 200);
con13 = (3*(185*x(1,2)+145*x(2,2)+105*x(3,2)+65*x(4,2)) + y(2,2) - w(2,2) >= 240);
con14 = (2.4*(185*x(1,2)+145*x(2,2)+105*x(3,2)+65*x(4,2)) + y(3,2) - w(3,2) >= 240);
con15 = (w(2,3) + w(2,4) <= 20*(185*x(1,3)+145*x(2,3)+105*x(3,3)+65*x(4,3)));
con16 = (w(3,3) + w(3,4) <= 16*(185*x(1,3)+145*x(2,3)+105*x(3,3)+65*x(4,3)));
con17 = (w(2,3) <= 6000);
con18 = (w(3,3) <= 6000);

constraints = [con1; con2; con3; con4; con5; con6; con7; con8; con9; con10; con11; con12; con13; con14; con15; con16; con17; con18];

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

% compute expected value of perfect information
value0 = value(obj);
x = binvar(4,3);
y = sdpvar(1,2);
w = sdpvar(1,4);

% objective
obj = 170*w(1) + 150 *w(2) + 36*w(3) + 10*w(4) - 170*1.4*y(1) - 150*1.4*y(2) -... 
         150*(185*x(1,1)+145*x(2,1)+105*x(3,1)+65*x(4,1)) -...
         230 *(185*x(1,2)+145*x(2,2)+105*x(3,2)+65*x(4,2)) -...
         260*(185*x(1,3)+145*x(2,3)+105*x(3,3)+65*x(4,3));

% constraints
con2 = (2.5*(185*x(1,1)+145*x(2,1)+105*x(3,1)+65*x(4,1)) + y(1) - w(1) >= 200);
con3 = (3*(185*x(1,2)+145*x(2,2)+105*x(3,2)+65*x(4,2))+ y(2) - w(2) >= 240);
con4 = (w(3) + w(4) <= 20*(185*x(1,3)+145*x(2,3)+105*x(3,3)+65*x(4,3)));
con5 = (w(3) <= 6000);
con7 = (y >= 0 );
con8 = (w >= 0);
con9 = (sum(x(1,:)) == 1);
con10 = (sum(x(2,:)) == 1);
con11 = (sum(x(3,:)) == 1);
con12 = (sum(x(4,:)) == 1);
constraints = [con2; con3; con4; con5; con7; con8; con9; con10; con11; con12];

% solve
optimize(constraints, -obj);
value1 = value(obj);

% constraints
con2 = (2.5*1.2*(185*x(1,1)+145*x(2,1)+105*x(3,1)+65*x(4,1)) + y(1) - w(1) >= 200);
con3 = (3*1.2*(185*x(1,2)+145*x(2,2)+105*x(3,2)+65*x(4,2))+ y(2) - w(2) >= 240);
con4 = (w(3) + w(4) <= 20*1.2*(185*x(1,3)+145*x(2,3)+105*x(3,3)+65*x(4,3)));
con5 = (w(3) <= 6000);
con7 = (y >= 0 );
con8 = (w >= 0);
con9 = (sum(x(1,:)) == 1);
con10 = (sum(x(2,:)) == 1);
con11 = (sum(x(3,:)) == 1);
con12 = (sum(x(4,:)) == 1);
constraints = [con2; con3; con4; con5; con7; con8; con9; con10; con11; con12];

% solve
optimize(constraints, -obj);
value2 = value(obj);

% constraints
con2 = (2.5*0.8*(185*x(1,1)+145*x(2,1)+105*x(3,1)+65*x(4,1)) + y(1) - w(1) >= 200);
con3 = (3*0.8*(185*x(1,2)+145*x(2,2)+105*x(3,2)+65*x(4,2))+ y(2) - w(2) >= 240);
con4 = (w(3) + w(4) <= 20*0.8*(185*x(1,3)+145*x(2,3)+105*x(3,3)+65*x(4,3)));
con5 = (w(3) <= 6000);
con7 = (y >= 0 );
con8 = (w >= 0);
con9 = (sum(x(1,:)) == 1);
con10 = (sum(x(2,:)) == 1);
con11 = (sum(x(3,:)) == 1);
con12 = (sum(x(4,:)) == 1);
constraints = [con2; con3; con4; con5; con7; con8; con9; con10; con11; con12];

% solve
optimize(constraints, -obj);
value3 = value(obj);

ExpeValuePerfectInfor= mean([value1, value2, value3]) - value0

end