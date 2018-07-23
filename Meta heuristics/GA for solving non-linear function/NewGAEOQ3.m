function [xm,fv,saveinventory]=NewGAEOQ3  %没有约束条件，有供应商，基于k.Deeb的论文
%%初始目标函数与约束条件
clear;
D=500;p=1/1000;Dcoe=5;Ab=50;Av=80;hb=10;hv=6;backcost=200;a=1;b=4;c=1;  %%问题模型中的参数
a0=0;breal=0.15;bint=0.35;preal=10;pint=4; %%遗传算法中的参数
syms z k; 
intloss=matlabFunction(int((z-k)*normpdf(z,0,1),z,k,inf)); %定义缺货积分函数
f=@(q,k,m)(Av*D/(m*q)+hv*q*(m*(1-D*p)-1+2*D*p)/2+Ab*D/q+hb*(q/2+k*Dcoe*sqrt(a*(b-q)^2+c))+backcost*D*Dcoe*sqrt(a*(b-q)^2+c)*intloss(k)/q);
f=@(x)f(x(1),x(2),x(3));
qmin=1e-1;
qmax=10000;  %变量上下界
kmin=1e-1;
kmax=100;
mmin=1e-1;
mmax=100;
NP=100; %进化多少代，人工指定的值

%%初始化种群
size=30; %种群个数,种群个数为求解变量的个数乘以10
E=zeros(size,6); %前三列为初始解，第四列为适应函数值，第五列记录是否为可行解，第六列记录违背约束条件的差值
E(:,1)=qmin+(qmax-qmin)*rand(size,1);
E(:,2)=kmin+(kmax-kmin)*rand(size,1);
E(:,3)=ceil(mmin+(mmax-mmin)*rand(size,1));
fv=inf;%初始最优值为无穷大的值
D=zeros(NP,4);%用来记录每代的最优解,平均值，最差解,最优解是否为可行解

%%计算适应函数罚函数值
for i=1:size
    E(i,4)=f(E(i,1:3));
end
fmax=max(E(:,3));

%%遗传进化   %%到这步适应值还没出错
for g=1:NP    %%原来错误在这里，这个k跟前面的k重复了
    %for i=1:size   %%小生态技术
        %A=randperm(size,1);
        %dij=1;
        %j=0;
        %while j<=floor(0.25*size);
            %if dij<0.1 
               %B=randperm(size,1);
               %if E(A,r)<=E(B,4)
                   %M(i,:)=E(A,(1:3));
               %else
                   %M(i,:)=E(B,(1:3));
               %end
               %break;
            %else
                %B=randperm(size,1);
                %dij=sqrt(1/3*(E(A,1)-E(B,1)).^2/(qmax-qmin).^2+1/3*(E(A,2)-E(B,2)).^2/(kmax-kmin).^2+1/3*(E(A,3)-E(B,3)).^2/(mmax-mmin).^2);
                %j=j+1;
            %end
        %end
        %if j>floor(0.25*size)
           %M(i,:)=E(A,(1:3));
        %end
    %end
    
    M=zeros(size,3);%用来存储优胜者的中间矩阵    
    for i=1:size   %%竞标赛选择出父母代，这种竞标赛选择方式应该是对的，合理的
        A=randperm(size,2);  
        if E(A(1),4)<=E(A(2),4)
            M(i,:)=E(A(1),(1:3));
        else
            M(i,:)=E(A(2),(1:3));
        end 
    end

    %%拉普拉斯交换
    E=zeros(size,6);  %将种群数据清零，用来记录后代的数据
    for j=1:size/2
        A=randperm(size,2);
        r=rand();   %对第一个变量交叉
        u=rand();
        if r>0.5
            beita=a0+breal*log(u);
        else
            beita=a0-breal*log(u);
        end
        E(2*j-1,1)=M(A(1),1)+beita*abs(M(A(1),1)-M(A(2),1));
        E(2*j,1)=M(A(2),1)+beita*abs(M(A(1),1)-M(A(2),1));
    %%只在可行解时出错是怎么回事？是不是变异的原因，已纠正
        r=rand();   %对第二个变量交叉
        u=rand();
        if r>0.5
            beita=a0+breal*log(u);
        else
            beita=a0-breal*log(u);
        end
        E(2*j-1,2)=M(A(1),2)+beita*abs(M(A(1),2)-M(A(2),2));
        E(2*j,2)=M(A(2),2)+beita*abs(M(A(1),2)-M(A(2),2));
        r=rand();   %对第三个变量交叉
        u=rand();
        if r>0.5
            beita=a0+bint*log(u);
        else
            beita=a0-bint*log(u);
        end
        E(2*j-1,3)=M(A(1),3)+beita*abs(M(A(1),3)-M(A(2),3));
        E(2*j,3)=M(A(2),3)+beita*abs(M(A(1),3)-M(A(2),3));
    end

    %%鲍威尔变异，该变异应该不会导致可行解不可行，对Deb论文的一个改进
    for i=1:size
        s1=rand();
        s=s1^preal;
        t1=(M(i,1)-qmin)/(qmax-M(i,1)); %对第一个变量变异
        if t1<rand()  
            E(i,1)=E(i,1)-s*(E(i,1)-qmin);
        else
            E(i,1)=E(i,1)+s*(E(i,1)-qmin);
        end
        t2=(M(i,2)-kmin)/(kmax-M(i,2)); %对第二个变量变异
        if t2<rand()  
            E(i,2)=E(i,2)-s*(E(i,2)-kmin);
        else
            E(i,2)=E(i,2)+s*(E(i,2)-kmin);
        end
        t3=(M(i,3)-mmin)/(mmax-M(i,3)); %对第三个变量变异
        s=s1^pint;
        if t3<rand()  
            E(i,3)=E(i,3)-s*(E(i,3)-kmin);
        else
            E(i,3)=E(i,3)+s*(E(i,3)-kmin);
        end
    end
    %%对第三个变量取整操作
    for i=1:size
        if fix(E(i,3))~=E(i,3)
            v=rand();
            if v>0.5
                E(i,3)=fix(E(i,3))+1;
            else
                E(i,3)=fix(E(i,3));
            end
        end
    end
    
    %%计算子代罚函数值，判断是否满足可行解 
    for i=1:size
        E(i,4)=f(E(i,(1:3)));
    end
    %for i=1:size
        %B(1)=subs(g1,[q,k],E(i,(1:2)));
        %if B(1)>=0
            %E(i,4)=1;
            %E(i,3)=subs(f,[q,k],E(i,(1:2)));%%跟直接算的结果不一样，也跟EOQ得到的结果不一样 eval出错的原因
        %else
            %E(i,4)=0;
            %E(i,3)=0;
        %end
        %if B(1)>=0
            %B(1)=0;
        %else
            %B(1)=abs(B(1));
        %end
        %E(i,5)=B(1);
    %end
    %for i=1:size
        %if E(i,4)<1e-6
            %E(i,3)=fmax+E(i,5);
        %end
    %end       
    [Q,IX]=sort(E,1);
    %Q=vpa(Q,4);%%
    D(g,1)=Q(1,4); %记录每代的最优值
    D(g,2)=mean(Q(:,4));
    D(g,3)=Q(size,4); 
    %D(g,4)=E(IX(1,3),4);
    if Q(1,4)<fv 
        xm=E(IX(1,4),(1:3));
        fv=Q(1,4);
        %xm=[xm,E(IX(1,3),4)];
        saveinventory=vpa(E(IX(1,3),2)*Dcoe*sqrt(a*(b-E(IX(1,3),1))^2+c),6);% 计算安全库存;
    end
end

%画图
t=1;
for i=20:NP
    %if D(i,4)==1
        N(t,:)=D(i,:);
        t=t+1;
    %end
end
        
plot(N(:,1),'r*');
hold on
plot(N(:,2),'b+');
hold on
plot(N(:,3),'ms');
legend('最优值','平均值','最差值');
hold off
D
end
   


 