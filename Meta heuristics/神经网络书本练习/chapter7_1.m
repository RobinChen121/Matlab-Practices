%%清空环境变量
clc
clear
%%产生输入 输出数据
% 设置步长
interval=0.01;
x1=-1.5:interval:1.5;
x2=-1.5:interval:1.5;
% 按照函数先求得相应的函数值，作为网络的输出
F=20+x1.^2-10*cos(2*pi*x1)+x2.^2-10*cos(2*pi*x2);
% %网络建立和训练
% 网络建立 输入为[x1;x2]，输出为F. Spread使用默认
net=newrbe([x1;x2],F)