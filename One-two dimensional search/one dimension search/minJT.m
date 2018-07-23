function [low,up]=minJT(fun,x0,h)
% fun是一个匿名函数

lamda=1;
x1=x0;x2=x0+h;

if fun(x1)>fun(x2)
    x_panduan=x2;
    while 1
        x1=x2;x2=x1+lamda*h;
        lamda=1.1*lamda;
        if fun(x2)>fun(x_panduan)
            break;
        end
    end
else
    x_panduan=x1;
    while 1
        x1=x2;x2=x1-lamda*h;
        lamda=1.1*lamda;
        if fun(x2)>fun(x_panduan)
            break;
        end
    end
end

low=min(x_panduan,x2);
up=max(x_panduan,x2);

end