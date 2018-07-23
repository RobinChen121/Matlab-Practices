function [x,minf]=miFactor(f,x0,g,w,alpha,var,eps)
%目标函数：f;
%初始点：x0;
%约束函数：g;
%乘子：w;
%放大系数：alpha；
%自变量向量：var;
%精度：eps；
%目标函数取最小值时的自变量值：x;
%目标函数的最小值：minf

format long;
if nargin==6
    eps=1.0e-4;
end

n=length(g);
while l
    for i=1:n
        gv(i)=subs(g,var,x0);
        if gv(i)>w(i)/alpha
            fg(i)=-0.5*w(i)^2/alpha;
        else
            fg(i)=0.5*((w(i)-alpha*g(i))^2-w(i)^2)/alpha;
        end
    end
    newf=f+sum(fg);
    [xm,minf]=minTruA(newf,x0,0.1,0.3,0.7,var);%用信赖域方法求解
    

FE=0;
