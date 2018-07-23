function carbon_ga(CityNum,shu)
%%shu=1表示计算碳成本，shu=0表示不计算碳成本
    [dislist,Clist]=tsp1(CityNum);
    DD=demand(CityNum);
global C
global n
global pc
global pm
global inn
global g
global P
global miu
global cishu
global M
global G_M
global G_M0
global carbon_per;
global tax
global a %%机会成本参数

inn=100; %初始种群大小
gnmax=100;  %最大代数
if CityNum==10
    C=10;
else
	C=20; %每辆车的载重量
end
g=50000;%每辆车的固定运营成本
P=150; %碳配额购买价
M=1000;% 碳排放限额
G_M=500;%每年的固定碳排放量，单位吨
G_M0=350;
n=floor(sum(DD)/C)+1; %车辆数
miu=7.7; %油价 
cishu=50; %年运货次数
pc=0.8; %交叉概率
pm=0.02; %变异概率
carbon_per=3.0772;%碳排放系数  kg/
tax=200; %碳税价格
a=5000;
 
ymax=zeros(gnmax,1);%用来记录每代的最优值
ymean=zeros(gnmax,1);%用来记录每代的平均值
smax=zeros(gnmax,CityNum);%用来记录每代的最优解
kx_n=zeros(gnmax,1);%用来记录每代最优解是否为可行解
qs_n=zeros(gnmax,n+1);%用来记录每代最优解的车辆的出发点
load_n=zeros(gnmax,CityNum+1);%用来记录每代最优解的各路径载重
Tgc_n=zeros(gnmax,1);%用来记录每代最优解的碳排放量

%产生初始种群
s=zeros(inn,CityNum);
tt=zeros(inn,1); %记录适应度
kx=zeros(inn,1);%记录是不是可行解
over=zeros(inn,1);%记录非可行解的超出量
qs=zeros(inn,n+1);%用来记录每个车辆的出发点与终止点的
load=zeros(inn,CityNum+1);%用来记录各路径载重
Tgc=zeros(inn,1);

for i=1:inn
    s(i,:)=randperm(CityNum);
    [tt(i),kx(i),qs(i,:),over(i),load(i,:)]=tcost(s(i,:),dislist,DD,shu);  %得到每个种群的适应值，并记录是否可行解
end   %load出错了
%%修改非可行解的适应值，并计算累计概率
max=0; %用来寻找可行解的最小值
for i=2:inn
    if  max<tt(i)&&kx(i)==1
        max=tt(i);
    end
end
for i=1:inn
    if kx(i)==0
        tt(i)=max+over(i)*tt(i)/C;%%修改非可行解
    end
end
[BB,idx]=sort(tt,1);
 %%开始计算累计概率,用轮盘赌方法选出
sumb=sum(BB);
pp=zeros(inn,1);
pp(1)=BB(1)/sumb;
smnew=zeros(inn,CityNum);%用来记录新矩阵
 for i=1:inn-1
     pp(i+1)=BB(i+1)/sumb+pp(i); 
 end
%%遗传操作
gn=1;
ymax(gn)=BB(1);
smax(gn,:)=s(idx(1),:);

while gn<gnmax+1
    for i=1:2:inn
        seln1=sel(s,idx,pp);  %选择父亲
        seln2=sel(s,idx,pp);  %选择母亲
        scro=cro(s,seln1,seln2,pc);  %交叉操作
        smnew(i,:)=mut(scro(1,:),pm);  %变异操作
        smnew(i+1,:)=mut(scro(2,:),pm);
    end
   s=smnew;  %产生了新的种群
   %%计算新种群适应值
   for i=1:inn
    s(i,:)=randperm(CityNum);
    [tt(i),kx(i),qs(i,:),over(i),load(i,:),Tgc(i,:)]=tcost(s(i,:),dislist,DD,shu);  %得到每个种群的适应值，并记录是否可行解
   end
    %%修改非可行解的适应值，并计算累计概率
    max=0; %用来寻找可行解的最大值
    for i=1:inn
        if  max<tt(i)&&kx(i)==1
            max=tt(i);
        end
    end
    for i=1:inn
        if kx(i)==0
            tt(i)=max+over(i)*tt(i)/C;
        end
    end
    [BB,idx]=sort(tt,1);
    %%开始计算累计概率,用轮盘赌方法选出
    sumb=sum(BB);
    pp=zeros(inn,1);
    pp(1)=BB(1)/sumb;
    smnew=zeros(inn,CityNum);
    for i=1:inn-1
        pp(i+1)=BB(i+1)/sumb+pp(i); 
    end
    %记录当前代的最佳个体
    ymax(gn)=BB(1);
    ymean(gn)=mean(BB(:,1));
    smax(gn,:)=s(idx(1),:);
    qs_n(gn,:)=qs(idx(1),:);
    kx_n(gn,:)=kx(idx(1),:);
    load_n(gn,:)=load(idx(1),:);
    Tgc_n(gn,:)=Tgc(idx(1),:);
    drawTSP10(Clist,smax(gn,:),qs(gn,:),ymax(gn),gn,0);
    gn=gn+1;
end
[BB,idx]=sort(ymax,1);
drawTSP10(Clist,smax(idx(1),:),qs_n(idx(1),:),ymax(idx(1)),idx(1),1);
smaxx=smax(idx(1),:);
load_nn=load_n(idx(1),:);
disp('起始点：');
qs_nn=smaxx(qs_n(idx(1),1:n))
disp('运输中的碳排放量：');
Tgc_nn=Tgc_n(idx(1),:)
disp('碳交易成本：');
P*(Tgc_nn-M)
disp('碳税成本：');
tax*Tgc_nn
disp('运输总成本');
BB(1)
disp('碳税总成本');
BB(1)-P*(Tgc_nn-M)+tax*Tgc_nn


figure(2);
plot(ymax,'r'); hold on;
plot(ymean,'b');grid; %%grid是为了显示网络线
title('搜索过程');
legend('最优解','平均解');
end

function [tt1,kx,qs,over,load,Tgc]=tcost(s,dislist,DD,shu)  %计算适应值
    m=1; %记录所用车辆数
    global n
    global C
    global miu
    global g
    global cishu
    global M
    global P
    global tax
    global carbon_per
    global G_M
    global G_M0
    global a
    qs=zeros(1,n+1);%记录车辆起始点 
    gf=zeros(1,n);%记录每辆车的总油耗成本
    gc=zeros(1,n);%记录每辆车的碳排放量
    CityNum=size(DD,2);
    sum1=0;
    qs(m)=1;%第一辆车的起始点为第一个点
    qs(n+1)=CityNum+1;%最后一辆车行驶的最后一个点
    for j=1:CityNum
        sum1=DD(s(j))+sum1;
        if sum1>C && m<n   %%判断可行解有点问题,没问题
            sum1=DD(s(j));
            m=m+1;
            qs(m)=j;
        end
    end
    load=zeros(1,CityNum+1);%计算行驶在每个路段上车辆的载重量
    m=1;
    for i=1:CityNum
        load(i)=sum(DD(s(qs(m):qs(m+1)-1)))-sum(DD(s(qs(m):i-1))); %没有问题
       if i+1>=qs(m+1)&&m<n
            m=m+1;
        end
    end
    m=1;
    while m<n+1
          sum11=miu*dislist(s(qs(m+1)-1),CityNum+1)*fuel_quantity(0)+miu*dislist(CityNum+1,qs(m))*fuel_quantity(load(qs(m)));       
          sum12=carbon_per*dislist(s(qs(m+1)-1),CityNum+1)*fuel_quantity(0)+2.6754*dislist(CityNum+1,qs(m))*fuel_quantity(load(qs(m)));
          for j=qs(m):qs(m+1)-2
              sum11=miu*dislist(s(j),s(j+1))*fuel_quantity(load(j+1))+sum11;%记录单辆车所走路段的总成本
              sum12=carbon_per*dislist(s(j),s(j+1))*fuel_quantity(load(j+1))+sum12;%记录单辆车所走路段的总碳排放量    
          end
          gf(m)=sum11;
          gc(m)=sum12/1000; %将单位化为吨
          m=m+1;
    end
    Tgc=cishu*sum(gc); % sum(gc)为运输总碳排放量
    tt1=cishu*sum(gf)*miu+n*g;
    if shu==1
        tt2=P*max((cishu*sum(gc)+G_M-M),0);
    end
    if shu==0
        tt2=0;
    end
    if shu==1
        tt2=tax*max((cishu*sum(gc)+G_M),0);
    end
    if shu==0
        tt2=0;
    end
    tt=tt1+tt2;
    if load(qs(m-1))<=C  %%有问题，少个等号
        kx=1; 
        over=0;
    else
        kx=0;
        over=load(qs(m-1))-C;
    end
end

function ff=fuel_quantity(load) %计算载重量的耗油
    global C
    if C==20
        ff=10.2181+1.021*load-1.1247*10^(-5)*load^2;
    else
        ff=6.5526+1.2154*load-3.5404*10^(-6)*load^2;
    end
end

function seln=sel(s,idx,pp) %选择
ppr=rand;
global inn
    for j=1:inn
        if pp(j)>=ppr
            seln(1,:)=s(idx(j),:);  %选择父亲
            break
        end
    end
end

function scro=cro(s,seln1,seln2,pc)
bn=size(s,2); %城市数
scro(1,:)=seln1;
scro(2,:)=seln2;%用来记录交叉后的路径
if rand<=pc
   c1=round(rand*(bn-2))+1;  %在[1,bn-1]范围内随机产生一个交叉位
   c2=round(rand*(bn-2))+1;
   chb1=min(c1,c2);
   chb2=max(c1,c2);
   middle=scro(1,chb1+1:chb2);
   scro(1,chb1+1:chb2)=scro(2,chb1+1:chb2);
   scro(2,chb1+1:chb2)=middle;
   for i=1:chb1       %找重复的点，替换掉
       for j=chb1+1:chb2
           if scro(1,j)==scro(1,i)
              scro(1,i)=scro(2,j);
           end
           if scro(2,j)==scro(2,i)
              scro(2,i)=scro(1,j);
           end
       end
   end
   for i=chb2+1:bn
       for j=1:chb2
           if scro(1,j)==scro(1,i)
              scro(1,i)=scro(2,j);
           end
           if scro(2,j)==scro(2,i)
              scro(2,i)=scro(1,j);
           end
       end
   end
end
end

%“变异”操作
function snnew=mut(snew,pm)

bn=size(snew,2);
snnew=snew;
if rand<pm
   c1=round(rand*(bn-2))+1;  %在[1,bn-1]范围内随机产生一个变异位
   c2=round(rand*(bn-2))+1;
   chb1=min(c1,c2);
   chb2=max(c1,c2);
   x=snew(chb1+1:chb2);
   snnew(chb1+1:chb2)=fliplr(x);%将这部分倒置
end
end

