function [minx,maxx]=minJT(f,x0,h0,eps)
%黄金分割法
format long;
if nargin==3 % nargin是用来判断输入变量个数的，若为3个，则默认误差精度
    eps=1.0e-6;
end

x1=x0;
k=0;
h=h0;
while 1
    x4=x1+h;
    k=k+1;
    f4=subs(f,findsym(f),x4);
    f1=subs(f,findsym(f),x1);
    if f4<f1
        x2=x1;
        x1=x4;
%         f2=f1;
%         f1=f4;
        h=2*h;
    else
        if k==1
            h=-h;
            x2=x4;
%             f2=f4;
        else
            x3=x2;
            x2=x1;
            x1=x4;
            break;
        end
    end
end
minx=min(x1,x3);
maxx=x1+x3-minx;
format short;