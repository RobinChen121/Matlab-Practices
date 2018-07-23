function [x,minf]=minZoujk(f,g,x0,var,eps)
% f为目标函数
% g为约束条件
% x0为初始点
% var为求解变量
% eps为精度要求

format long;
if nargin==4;
    eps=1.0e-6;
end

syms lamda;
n=length(g);
m=length(var);
df=jacobian(f,var);
c=zeros(1,m+1);
c(m+1)=1;
A=zeros(n+1,m+1);
b=zeros(n+1,1);
lb=zeros(m+1,1);
ub=zeros(m+1,1);
j=1;
for i=1:m
    lb(i)=-1;
    ub(i)=1;
end
lb(m+1)=-inf;
ub(m+1)=inf;

for i=1:n
    dg(i,:)=jacobian(g(i),var);
end

while 1
    fv=subs(df,var,x0);
    for i=1:n
        dgv(i,:)=subs(dg(i,:),var,x0);
        gv(i,:)=subs(g(i,:),var,x0);
    end
    b(1)=0;
    for i=2:(n+1)
        b(i)=gv(i-1);
    end
    A(1,:)=[fv,-1];
    for i=2:(n+1)
        A(i,:)=[-dgv(i-1,:),-1];
    end
    zk=linprog(c,A,b,[],[],lb,ub);
    if abs(zk(m+1))<=1.0e-6
        break;
    else
        dk=transpose(zk(1:m));
        x1=x0+lamda*dk;
        tmpf=subs(f,var,x1);
        for i=1:n
            gv1(i)=subs(g(i,:),var,x1);
            mlamda(i)=solve(gv1(i),lamda);
        end
        maxlamda=max(double(mlamda));    
        [xm,minf]=minHJ(tmpf,0,maxlamda,1.0e-14);%用黄金分割法求解一维极值问题
        x0=x0+xm*dk;
    end
    Tg(j)=subs(f,var,x0)
    j=j+1;
end
x=x0;
minf=subs(f,var,x0);
format short;
        
    
    
    
    