function [xm,fv]=GAEOQ1
%%初始目标函数与约束条件
%求解变量t, z, 目标函数 f,约束条件 g1
syms t k z;
f=100/t+25*(25*t+k*10*sqrt(1+1.25*t))+100*10*sqrt(1+1.25*t)*int((z-k)*normpdf(z,0,1),z,k,inf)/t;
tmin=1e-2;
tmax=4;
kmin=0;
kmax=2;
g1=1-200*sqrt(2+0.75*t)*int((z-k)*normpdf(z,0,1),z,k,inf)/(50*t);%%将约束标准化，将右端变为1
NP=50; %进化多少代

%%初始化种群，种群长度40
size=20;
E=zeros(20,5); %前两列为初始解，第三列为适应函数值，第四列记录是否为可行解，第五列记录违背约束条件的差值
E(:,1)=tmin+(tmax-tmin)*rand(size,1);
E(:,2)=kmin+(kmax-kmin)*rand(size,1);
fv=inf;%初始最优值为无穷大的值
D=zeros(NP,4);%用来记录每代的最优解,平均值，最差解,最优解是否为可行解


%%计算适应函数罚函数值，判断是否为可行解
for i=1:size
    B=zeros(1,1);
    B(1)=subs(g1,[t,k],E(i,(1:2)));
    if B(1)>=0
       E(i,4)=1;
       E(i,3)=subs(f,[t,k],E(i,(1:2)));
    else
       E(i,4)=0;
       E(i,3)=0;
    end
    if B(1)>=0
       B(1)=0;
    else
       B(1)=abs(B(1));
    end
    E(i,5)=B(1);
end
fmax=max(E(:,3));
for i=1:size
if E(i,4)<1e-6
   E(i,3)=fmax+E(i,5);
   end
end

%%遗传进化   %%到这步适应值还没出错
for g=1:NP    %%原来错误在这里，这个k跟前面的k重复了
    %%竞标赛选择  %%小生态技术
    M=zeros(size,2);%用来存储优胜者的中间矩阵
        for i=1:size
            %A=randperm(size,6);
            %dij1=sqrt(0.5*(E(A(1),1)-E(A(2),1))^2);  %%小生态技术,只有单界时小生态技术没法用
            %dij2=sqrt(0.5*(E(A(1),1)-E(A(3),1))^2); 
            %dij3=sqrt(0.5*(E(A(1),1)-E(A(4),1))^2);
            %dij4=sqrt(0.5*(E(A(1),1)-E(A(5),1))^2);
            %dij5=sqrt(0.5*(E(A(1),1)-E(A(6),1))^2);
            %if dij1<0.1
                %if E(A(1),3)<=E(A(2),3)
                    %M(i,:)=E(A(1),(1:2));
                %else
                    %M(i,:)=E(A(2),(1:2));   
                %end
                %continue;
            %elseif dij2<0.1
                %if E(A(1),3)<=E(A(3),3)
                    %M(i,:)=E(A(1),(1:2));
                %else
                    %M(i,:)=E(A(3),(1:2));   
                %end
                %continue;
            %elseif dij3<0.1
                %if E(A(1),3)<=E(A(4),3)
                    %M(i,:)=E(A(1),(1:2));
                %else
                    %M(i,:)=E(A(4),(1:2));   
                %end
                %continue;
            %elseif dij4<0.1
                %if E(A(1),3)<=E(A(5),3)
                    %M(i,:)=E(A(1),(1:2));
                %else
                    %M(i,:)=E(A(5),(1:2));   
                %end
                %continue;
            %elseif dij5<0.1
                %if E(A(1),3)<=E(A(6),3)
                    %M(i,:)=E(A(1),(1:2));
                %else
                    %M(i,:)=E(A(6),(1:2));   
                %end
            %else
                %M(i,:)=E(A(1),(1:2));
            %end  
        %end
            A=randperm(size,2);
            if E(A(1),3)<=E(A(2),3)
                M(i,:)=E(A(1),(1:2));
            else
                M(i,:)=E(A(2),(1:2));
            end 
        end

    %%模拟二进制交叉生成后代
        for j=1:size/2
            if rand()>=0.5
                A=randperm(size,2);
                c=rand();
                x2=max(M(A(1),1),M(A(2),1));
                x1=min(M(A(1),1),M(A(2),1));
                beita1_t=1+2*(x1-tmin)/(x2-x1);
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
        end

        %%变异，变异会不会导致可行解不可行?单下界时不用变异,设置判断条件防止过界
        for i=1:size
            nita=100+g;
            pm=1/size+g*(1-1/size)/NP;
            if rand()<pm
                u=rand();
                %x=E(i,1);
            	x=E(i,1);%%
                deltamax=1;
                if u<=0.5
                    delta_2=(2*u)^(1/(nita+1))-1;
                else
                    delta_2=1-(2*(1-u))^(1/(nita+1));
                end
                if x+delta_2*deltamax>=tmin
                    E(i,1)=x+delta_2*deltamax;
                end
            end
    
            if rand()<pm
                u=rand();
                x=E(i,2);
                deltamax=1;
                if u<=0.5
                    delta_1=(2*u)^(1/(nita+1))-1;
                else
                    delta_1=1-(2*(1-u))^(1/(nita+1));
                end
                if x+delta_1*deltamax>=kmin
                    E(i,2)=x+delta_1*deltamax;
                end
            end
        end
    %%计算子代罚函数值，判断是否满足可行解 
    for i=1:size
        B(1)=subs(g1,[t,k],E(i,(1:2)));
        if B(1)>=0
            E(i,4)=1;
            E(i,3)=subs(f,[t,k],E(i,(1:2)));%%跟直接算的结果不一样，也跟EOQ得到的结果不一样 eval出错的原因
        else
            E(i,4)=0;
            E(i,3)=0;
        end
        if B(1)>=0
            B(1)=0;
        else
            B(1)=abs(B(1));
        end
        E(i,5)=B(1);
    end
    for i=1:size
        if E(i,4)<1e-6
            E(i,3)=fmax+E(i,5);
        end
    end       
    [Q,IX]=sort(E,1);
    %Q=vpa(Q,4);%%
    D(g,1)=Q(1,3);
    D(g,2)=mean(Q(:,3));
    D(g,3)=Q(size,3); 
    D(g,4)=E(IX(1,3),4);
    if Q(1,3)<fv && D(g,4)==1
        fv=Q(1,3);
        xm=E(IX(1,3),(1:2));
        xm=[xm,E(IX(1,3),4)];
    end
end

%画图
k=1;
for i=1:NP
    if D(i,4)==1
        N(k,:)=D(i,:);
        k=k+1;
    end
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
   


 