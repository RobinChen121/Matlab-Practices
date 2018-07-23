% myfun=@(x1,x2)(1+(x1+x2+1)^2*(19-14*x1^2+3*x1^2-14*x2+6*x1*x2+3*x2^2))*(30+(2*x1-3*x2)^2*(18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2))

% This is a standalone version of glbSolve.m which is a part of the
% optimization environment TOMLAB, see www.ima.mdh.se/tom

%gblSolve implepents the algorithm DIRECT by D.R.Jones, C.D. Perttune and
%B.E.Stuckman presented in the paper "Lipschitzian OPtimization Without the
%Lipschitz Constant", Vol.79, No.1, October 1993

%gblSolve solves problems of the form:
%
% min f(x)
% x
% s.t. x_L<=x<=x_U
%
%Calling syntax:
%function Result=gblSolve(fun,x_L,x_U,GLOBAL,Prilev)
%
%INPUT PARAMETERS
%
% fun  Name of m-file computing the function value, given as a string.
% x_L  Lower bounds for x
% x_U  Upper bounds ofr x
%
% GLOBAL.iterations    Number of iterations to run, default 50.
% GLOBAL.epsilon       Globla/local search weight parameter, default 1e-4.
%
% If restart is wanted, the following fields in GLOBAL should be defined
% and equal the corresponding fields in the Result structure from the
% previous run:
% GLOBAL.C   Matrix with all rectangle centerpoints.
% GLOBAL.D   Vector with distances from centerpoints to the vertices.
% GLOBAL.L   Matrix with all Rectangle side lengths in each dimension
% GLOBAL.F   Vector with function values.
% GLOBAL.d   Row vector of all different distances, sorted.
% GLOBAL.d_min   Row vector of minimum function value for each distance

% Prilev     Printing level:
%            PriLev>=0   Warnings
%            PriLev> 0   Small info
%            PriLev> 1   Each iteration info
%
% OUTPUT PARAMETERS
%
% Result     Structure with fields:
%    x_k     Matrix with all points fulfilling f(x)=min(f)
%    f_k     Smallest function value found.
%    Iter    Number of iterations
%    FuncEV  Number of function evaluations.
%    GLOBAL.C   Matrix with all rectangle centerpoints.
%    GLOBAL.D   Vector with distances from centerpoint to the vertices.
%    GLOBAL.L   Matrix with all rectangle side lengths in each dimension.
%    GLOBAL.F   Vector with function values.
%    GLOBAL.d   Row vector of all different distances, sorted.
%    GLOBAL.d_min   Row vector of minimum function value for each distance
%
%    written by Mattias Bjorkman, last modified May 4, 1999

function Result=gblSolve(fun,x_L,x_U,GLOBAL,PriLev)
if nargin<5
    PriLev=[];
    if nargin<4
        GLOBAL=[];
        if nargin<3
            x_U=[];
            if nargin<2
                x_L=[];
                if nargin<1
                    fun=[];
                end
            end
        end
    end
end

if isempty(PriLev)
    PriLev=1;
end
if isempty(fun)||isempty(x_L)||isempty(x_U)
    disp('gblSolve requires at least three nonempty input arguments');
    return;
end
if isempty(GLOBAL)
    T=50;   % Number of iterations
    epsilon=1e-4;   % global/local wight parameter.
    tol=0.01;       % Error tolerance parameter.
else
    if isfield(GlOBAL,'iterations')  % Number of iterations
        T=GLOBAL.iterations;
    else
        T=50;
    end
    if isfield(GLOBAL,'epsilon')  % global/local wight parameter
        epsilon=GLOBAL.epsion;
    else
        epsilon=1e-4;
    end
    if isfield(GLOBAL,'tolerance')  %Convergence tolerance
        tol=GLOBAL.tolerance;
    else
        tol=0.01;
    end
end

nFunc=0;
convflag=0;

x_L=x_L(:);
x_U=x_U(:);
n=length(x_L);  % Problem dimension 问题维度

tolle=1e-16;
tolle2=1e-12;

%
% step 1, Initializaion
%
if isfield(GLOBAL,'C')&&~isempty(GLOBAL.C)
    % Restart with values from previous run.
    F=GLOBAL.F;
    m=length(F);
    if PriLev>0
        fprintf('\n Restarting with %d sampled points from previous run\n',m);
    end
    
    D=GLOBAL.D;
    L=GLOBAL.L;
    d=GlOBAL.d;
    d_min=GLOBAL.d_min;
    
    f_min=min(F);
    E=max(epsilon*abs(f_min),1e-8);
    [dummy i_min]=min((F-f_min+E)./D);
    
    % Must transform GLOBAL.C back to unit hypercube
    for i=1:m
        C(:,i)=(GLOBAL.C(:,i)-x_L)./(x_U-x_L);
    end
else
    % No restart, set first point to center of the unit hypercube.
    m=1;  % Current number of rectangles
    C=ones(n,1)./2;  % Matrix with all rectangle centerpoints
    % All C_coordinates refers to the n-dimensional hypercubes.
    
    x_m=x_L+C.*(x_U-x_L); % Transform C to original search space  中点的坐标
    f_min=feval(fun,x_m); % Function value at x_m
    f_0=f_min;
    nFunc=nFunc+1;        % nFunc 用来统计 计算函数值 的次数
    i_min=1; % The rectangle which minimizes (F-f_min+E)./D where
             % E=max(epsilon*abs(f_min),1E-8)
    
    L=ones(n,1)./2;     % Matrix with all rectangle side lengths in each dimension
    D=sqrt(sum(L.^2));  % Vector with distances from centerpoint to the vertices
    F=[f_min];          % Vector with function values
    
    d=D;                % Row vector of all different distances,sorted
    d_min=f_min;        % Row vector of minimum function value for each distance
end

% ITERATION LOOP
t=1;    % t is the iteration counter
while t<=T && ~ convflag
    
    %
    % step 2  Identify the set S of all potentially optimal rectangles
    %
    S=[];    % Set of all potentially optimal rectangles
    
    S_1=[];
    idx=find(d==D(i_min));
    % 或者写成 idx=find(abs(d-D(i_min))<=tolle);
    if isempty(idx)
        if PriLev>=0
            fprint('\n Warning: Numerical trouble when determining S_1\n');
        end
        return;
    end
    for i=idx:length(d)
        idx2=find( ( F==d_min(i) )&( D>=d(i) ) ); % 纵坐标、横坐标
        % 或 idx2=find((abs(F-d_min(i))<=tolle)&& (abs(D-d(i))<=tolle));
        S_1=[S_1 idx2];   % S_1为从最小纵坐标到最大纵坐标的 index
    end
    % S_1 now includes all rectangles i, with D(i)>=D(i_min)
    % S_1包括了所有距离大于等于最小距离的点
    % and F(i) is the minimum function value for the current distance.
    
    % Pick out all rectangles in S_1 which lies below the line passing through
    % the points: (D(i_min),F(i_min)) and the lower rightmost point.
    S_2=[];
    if length(d)-idx>1
        a1=D(i_min);
        b1=F(i_min);
        a2=d(length(d));        % d是中点到顶点距离，a2是横坐标最大的
        b2=d_min(length(d));    % d_min是中点的函数值
        % The line is defined by: y= slope*x +const
        slope=(b2-b1)/(a2-a1);
        const=b1-slope*a1;      % 第一条线
        for i=1:length(S_1)
            j=S_1(i);
            if F(j)<=slope*D(j)+const+tolle2
                S_2=[S_2 j];
            end
        end
        % S_2 now contains all points in S_1 which les on or below the line
        
        % Find the points on the convex hull defined by the points in S_2
        xx=D(S_2);
        yy=F(S_2);
        h=conhull(xx,yy);  % conhull is an intern subfunction
        S_3=S_2(h);
    else
        S_3=S_1;
    end
    S=S_3;
    
    
    % step 3-5: Select any rectangle j in S
    for jj=1:length(S)  % For each potentially optimal rectangle
        j=S(jj);
        
        %
        % Step 4 determin where to sample within rectangle j and how to
        % devide the rectangle into subrectangles. Update f_min and set
        % m=m+delta_m, where delta_m is the number of new points sampled.
        
        % 4:1 Identify the set I of dimensions with the maximum side
        % length.
        %   let delta equal one-third of this maximum side length.
        
        max_L=max(L(:,j));
        I=find(L(:,j)==max_L);   % I　为最大边长的矩形
        % I=find(abs(L(L,j)-max_L)<tolle);
        delta=2*max_L/3;
        
        % 4:2 Sample the function at the points c +- delta*e_i for all i in
        % I
        w=[];
        for ii=1:length(I) % for each dimension with maximum side length
            i=I(ii);
            e_i=[zeros(i-1,1);1;zeros(n-i,1)];
            
            c_m1=C(:,j)+delta*e_i;       % Centerpoint for new rectangle
            x_m2=x_L+c_m2.*(x_U-x_L);    % Transform c_m2 to original search space
            f_m2=feval(fun,x_m2);        % Function value at x_m2
            nFunc=nFunc+1;
            
            c_m2=C(:,j)-delta*e_i;       % Centerpoint for new rectangle
            x_m2=x_L+c_m2.*(x_U-x_L);    % Transform c_m2 to original search space
            f_m2=feval(fun,x_m2);        % Function value at x_m2
            nFunc=nFunc +1;
            
            w(ii)=min(f_m1,f_m2);
            
            C=[C c_m1 c_m2];             % Matrix with all rectangle centerpoints
            F=[F f_m1 f_m2];             % Vector with function values
        end
        
        %  4:3 Divide the rectangle containing C(:,j) into thirds along the
        %  dimension in I, starting with dimension with the lowest value of
        %  w(ii)
        [a,b]=sort(w);
        for ii=1:length(I)
            i=I(b(ii));
            
            ix1=m+2*b(ii)-1;  % index for new rectangle
            ix2=m+2*b(ii);    % index for new rectangle
            
            L(i,j)=delta/2;
            
            L(:,ix1)=L(:,j);
            L(:,ix2)=L(:,j);
            
            D(j)=sqrt(sum(L(:,j).^2));
            D(ix1)=D(j);
            D(ix2)=D(j);
        end
        m=m+2*length(I);
    end
    
    % Update
    f_min=min(F);
    E=max(epsilon*abs(f_min),1e-8);
    [dummy i_min]=min((F-f_min+E)./D);
    
    d=D;
    i=1;
    while 1
        d_temp=d(i);
        idx=find(d~=d_tmp);
        d=[d_tmp d(idx)];
        if i==length(d)
            break;
        else
            i=i+1;
        end
    end
    d=sort(d);
        
    d_min=[];
    for i=1:length(d)
        idx1=find(D==d(i));
        % idx1=find(abs(D-d(i))<=tolle);
        d_min(i)=min(F(idx1));
    end
    
    if PriLev>1
        fprintf('\n Iteration: %d f_min: %15.10f  Sampled points: %d',t,f_min,nFunc);
    end
    
    t=t+1;
end  % Iteration loop;


% save results

Result.f_k=f_min;   % best function value
Result.Iter=T;      % Number of iterations

CC=[];
for i=1:m           % Transform to original coordinates
    CC=[CC x_L+C(:,i).*(x_U-x_L)];
end

Result.GlOBAL.C=CC;      % all sampled points in original coordinates
Result.GLOBAL.F=F;       % all function values computed
Result.GLOBAL.D=D;       % all distances
Result.GLOBAL.L=L;       % all lengths
Result.GLOBAL.d=d;
Result.GLOBAL.d_min=d_min;

% Find all points i with F(i)=f_min
idx=find(F==f_min);
Result.x_k=CC(:,idx); % all points i with F(i)=f_min;

Result.FuncEv=nFunc;
end


function h=conhull(x,y)   
% conhull returns all points on the convex hull, even redundant ones.
%
% conhull is based on the algorithm GRAHAMSHULL pages 108-109
% in "Computational Geometry" by Franco P. Preparata and
% Michael Ian Shamos.
%
% Input vector x must be sorted i.e. x(1)<=x(2)<=...<=x(length(x)).
%
x=x(:);
y=y(:);
m=length(x);
if length(x)~=length(y)
    disp('Input dimension must agree, error in conhull-gblSolve');;
    return;  % 现在的matlab版本 return 可以不用，自动返回
end
if m==2
    h=[1 2];
end
if m==1
    h=1;
end
v=1;             % 初始点为第一个点
w=length(x);     % 凸包中的点的个数
flag=0;
h=[1:length(x)]'; % Index vector for points in convex hull
while (next(v,m)~=1) || (flag==0)
    if next(v,m)==w  % flag的作用主要为了在出现这种情况时迭代终止：下一个点为最后一个点
        flag=1;
    end
    a=v;b=next(v,m);c=next(b,m);
    if det([x(a) y(a) 1; x(b) y(b) 1; x(c) y(c) 1])>=0
        leftturn=1;
    else
        leftturn=0;
    end
    if leftturn
        v=next(v,m);
    else
        x=[x(1:b-1);x(b+1:m)];
        y=[y(1:b-1);y(b+1:m)];
        h=[h(1:b-1);h(b+1:m)];   %   将b点排除
        m=m-1;
        w=w-1;
        v=pred(v,m);             %   a点退后一位
    end
end
end

function i=next(v,m)
if v==m
    i=1;
else
    i=v+1;
end
end

function i=pred(v,m)
if v==1
    i=m;
else
    i=v-1;
end
end



    
        
    
    

