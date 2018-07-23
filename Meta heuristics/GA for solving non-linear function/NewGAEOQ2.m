function [xm,fv,saveinventory]=NewGAEOQ2  %没有约束条件，有供应商，基于k.Deb的论文
%%初始目标函数与约束条件
clear;
D=1000;sigma=5;P=3200;p=1/P;K=400;F=25;A=50;hb=5;hv=4;b=0.01;pi1=100;
syms z k;
intloss=matlabFunction(int((z-k)*normpdf(z,0,1),z,k,inf));
syms q n;
Hn=hb+hv*(n*(1-D*p)-1+2*D*p);
Gn=F+(A+K)/n;
f=Gn*D/q+q*Hn/2+hv*k*sigma*sqrt(p*q+b)+pi1*D*sigma*sqrt(p*q+b)*intloss(k)/q;
%f=@(q,k,m)(Av*D/(m*q)+hv*q*(m*(1-D*p)-1+2*D*p)/2+Ab*D/q+hb*(q/2+k*Dcoe*sqrt(a*(b-q)^2+c))+backcost*D*Dcoe*sqrt(a*(b-q)^2+c)*intloss(k)/q);

qmin=b;  %%待改
qmax=10000;  %变量上下界
kmin=1e-2;
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
    E(i,4)=subs(f,[q,k,n],E(i,(1:3)));
end
fmax=max(E(:,3));

%%遗传进化   %%到这步适应值还没出错
for g=1:NP    %%原来错误在这里，这个k跟前面的k重复了
    %%竞标赛选择  %%小生态技术
    M=zeros(size,3);%用来存储优胜者的中间矩阵
    %for i=1:size
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
        
    for i=1:size   %%竞标赛选择出父母代
        A=randperm(size,2);  
        if E(A(1),4)<=E(A(2),4)
            M(i,:)=E(A(1),(1:3));
        else
            M(i,:)=E(A(2),(1:3));
        end 
    end

    %%模拟二进制交叉生成后代
    for j=1:size/2
        if rand()>=0.5 
            A=randperm(size,2);
            c=rand();
            x2=max(M(A(1),1),M(A(2),1));
            x1=min(M(A(1),1),M(A(2),1));
            beita1_t=1+2*(x1-qmin)/(x2-x1);
            rfa_t=2-beita1_t^(-2);
            if c<=1/rfa_t
                beita2_t=sqrt(rfa_t*c);
            else
                beita2_t=sqrt(1/(2-rfa_t*c));
            end
            E(2*j-1,1)=0.5*(x1+x2-beita2_t*(x2-x1));
            E(2*j,1)=0.5*(x1+x2+beita2_t*(x2-x1));
        end
    %%只在可行解时出错是怎么回事？是不是变异的原因，已纠正
        if rand()>0.5
            c=rand();
            x2=max(M(A(1),2),M(A(2),2));
            x1=min(M(A(1),2),M(A(2),2));
            beita1_t=1+2*(x1-kmin)/(x2-x1);
            rfa_t=2-beita1_t^(-2);
            if c<=1/rfa_t
                beita2_t=sqrt(rfa_t*c);
            else
                beita2_t=sqrt(1/(2-rfa_t*c));
            end
            E(2*j-1,2)=0.5*(x1+x2-beita2_t*(x2-x1));
            E(2*j,2)=0.5*(x1+x2+beita2_t*(x2-x1));
        end
         if rand()>0.5
            c=rand();
            x2=max(M(A(1),3),M(A(2),3));
            x1=min(M(A(1),3),M(A(2),3));
            beita1_t=1+2*(x1-mmin)/(x2-x1);
            rfa_t=2-beita1_t^(-2);
            if c<=1/rfa_t
                beita2_t=sqrt(rfa_t*c);
            else
                beita2_t=sqrt(1/(2-rfa_t*c));
            end
            E(2*j-1,3)=ceil(0.5*(x1+x2-beita2_t*(x2-x1)));
            E(2*j,3)=ceil(0.5*(x1+x2+beita2_t*(x2-x1)));
         end
    end

    %%变异，变异会不会导致可行解不可行?单下界时不用变异,设置判断条件防止过界
    for i=1:size
        nita=100+g;
        pm=1/size+g*(1-1/size)/NP;
        if rand()<pm
            u=rand();
            x=E(i,1);%%
            deltamax=qmax-qmin;
            delta_1=min(x-qmin,qmax-x)/deltamax;
            if u<=0.5
                delta=(2*u+(1-2*u)*(1-delta_1)^(1+nita))^(1/(nita+1))-1;
            else
                delta=1-(2*(1-u)+2*(u-0.5)*(1-delta_1)^(1+nita))^(1/(nita+1));
            end
            if x+delta*deltamax>=qmin %这是为了防止出现因变异而超出上下界的情况
                E(i,1)=x+delta*deltamax;
            end  
            x=E(i,2);
            deltamax=kmax-kmin;
            delta_1=min(x-kmin,kmax-x)/deltamax;
            if u<=0.5
                delta=(2*u+(1-2*u)*(1-delta_1)^(1+nita))^(1/(nita+1))-1;
            else
                delta=1-(2*(1-u)+2*(u-0.5)*(1-delta_1)^(1+nita))^(1/(nita+1));
            end
            if x+delta*deltamax>=kmin
                E(i,2)=x+delta*deltamax;
            end
            x=E(i,3);
            deltamax=mmax-mmin;
            delta_1=min(x-mmin,mmax-x)/deltamax;
            if u<=0.5
                delta=(2*u+(1-2*u)*(1-delta_1)^(1+nita))^(1/(nita+1))-1;
            else
                delta=1-(2*(1-u)+2*(u-0.5)*(1-delta_1)^(1+nita))^(1/(nita+1));
            end
            if x+delta*deltamax>=mmin
                E(i,3)=ceil(x+delta*deltamax);
            end
        end
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
        E(i,4)=subs(f,[q,k,n],E(i,(1:3)));
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
    D(g,1)=Q(1,4);
    D(g,2)=mean(Q(:,4));
    D(g,3)=Q(size,4); 
    %D(g,4)=E(IX(1,3),4);
    if Q(1,4)<fv 
        xm=E(IX(1,4),(1:3));
        fv=Q(1,4);
        %xm=[xm,E(IX(1,3),4)];
        saveinventory=vpa(xm(2)*sigma*sqrt(p*xm(1)+b),6);
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
   


 