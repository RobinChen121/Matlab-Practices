function AllNeighborSearch(C)
%%% author: Chen Zhen
%%% 本程序从所有城市中选取第一个城市作为初始点
    if nargin<1,
       error('There are no input data!')
    end
    if ~isnumeric(C),
        error('The array C must be numeric!') 
    end
    if ~isreal(C),
        error('The array C must be real!') 
    end
    s=size(C); % size of array C
    if length(s)~=2,
        error('The array C must be 2D!') 
    end
    if s(1)<3,
        error('Must be not less than 3 cities!')
    end
tic; % 计算程序执行所用时间
CityNum= size(C,1);
DistanceList=zeros(CityNum,CityNum);
for i=1:CityNum     % 求解城市间的两两距离
    for j=1:CityNum
        DistanceList(i,j)=((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5;
    end
end

S0=[1,randperm(CityNum-1)+1]
Distance0=CalDistance(DistanceList,S0); % 计算初始路径的距离
S=zeros(1,CityNum);
S=S0;%初始解
Distance=Distance0;%初始距离
NeighborNum=ceil((CityNum-1)*(CityNum-2)/2);%邻居数目
Neighbor=zeros(NeighborNum,CityNum);%邻域
Swap=zeros(NeighborNum,2);%2-opt邻域交换对象
a=1;
lx=50; % 每个局部邻域的迭代步数

for i=1:CityNum     % 表示2-opt邻域所有的交换对象
    Swap(a:a+CityNum-i-2,1)=i+1;
    Swap(a:a+CityNum-i-2,2)=i+2:CityNum;
    a=a+CityNum-i-1;
end

step=1;% 迭代步数
while Distance<=Distance0
    PresentDistance(step)=CalDistance(DistanceList,S);%对当前路径距离存储
    PresentRoute(step,:)=S;
    for i=1:NeighborNum   % 对初始路径的节点用交换对象替代
        Neighbor(i,:)=S;
        Neighbor(i,Swap(i,1))=S(1,Swap(i,2));
        Neighbor(i,Swap(i,2))=S(1,Swap(i,1));
    end
    D=zeros(NeighborNum,2);
    for i=1:NeighborNum   % 求所有邻域中路径的距离
        D(i,1)=i;
        D(i,2)=CalDistance(DistanceList,Neighbor(i,:));
    end
    [B,IX]=sort(D(:,2)); % 对路径距离按从小到大排序
    if B(1)<Distance
        Distance=B(1);
        BetterDistance(step)=Distance;%对改进路径距离进行存储
        S=Neighbor(IX(1),:); 
        step=step+1;
    else
        BetterDistance(step)=Distance;
        Route=S
        ShortDistance=Distance
        break
    end
end
toc % 计算程序执行所用时间
for i=1:step
    DrawRoute(C,PresentRoute(i,:),PresentDistance(i),i); % 对路径变化过程画图
end
figure(2); % 对距离变化过程画图
plot(BetterDistance,'r'); hold on;
plot(PresentDistance,'b');
title('搜索过程');
legend('改进解','当前解');
hold off;
end
