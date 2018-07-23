function [x,minf]=minBFGS(f,x0,var,eps)
%目标函数：f;
%初始点：x0;
%自变量向量：var;
%目标函数取最小值时的自变量值：x;
%目标函数的最小值：minf

format long;
if nargin==3
    eps=1.0e-3;
end
x0=transpose(x0);%初始值的转置
n=length(var);%自变量个数
syms l;%l默认值为1
H=eye(n,n); %初始H矩阵   
gradf=jacobian(f,var);  %梯度矩阵
v0=Funval(gradf,var,x0);%梯度矩阵值
p=-H*transpose(v0);%初始下降方向
k=0;
i=0;

while 1
    k=0;
    v=Funval(gradf,var,x0); %将初始值代入梯度矩阵
    tol=norm(v) %梯度矩阵的行列式值
    i=i+1;
    Tg(i)=Funval(f,var,x0)
    if tol<=eps %计算终止条件
        x=x0;
        break;
    end
    y=x0+l*p;%迭代值
    yf=Funval(f,var,y);
    [a,b]=minJT(yf,0,0.1);%用进退法确定搜索区间
    xm=minGX(yf,a,b)%用黄金分割法确定搜索步长
    x1=x0+xm*p;%下一个迭代点
    vk=Funval(gradf,var,x1);%新点的梯度值
    tol=norm(vk);%精度判断标准
    if tol<=eps
        x=x1;
        break;
    end
    if k==n
        x0=x1;
        continue;
    else
        dx=x1-x0;
        dgf=vk-v;
        dgf=transpose(dgf);
        dxT=transpose(dgf);
        dgfT=transpose(dgf);
        mdx=dx*dxT;
        mdgf=dgf*dgfT;
        H=H+(1+dgfT*(H*dgf)/(dxT*dgf))*mdx/(dxT*dgf)-(dx*dgfT*H+H*dgf*dxT)/(dxT*dgf);
        p=-H*transpose(vk);
        k=k+1;
        x0=x1;
    end
    x0
end
minf=Funval(f,var,x);
i
plot(Tg)
format short;



