function [xm,fv]=GAEOQ3
%%初始目标函数与约束条件
syms t k1 k2 k3 z m1 m2 m3;
f=50/t+30/(m1*t)+20/(x2*t)+25/(x3*t)+25*(25*m1*t+k1*10*sqrt(1+1.25*t))+1000*sqrt(1+1.25*t)*int((z-k1)*normpdf(z,0,1),z,k1,inf)/(m1*t)+30*(30*m2*t+k2*12*sqrt(1+1.25*t))+120*12*sqrt(1+1.25*t)*int((z-k2)*normpdf(z,0,1),z,k2,inf)/(m2*t)+20*(20*m3*t+k3*15*sqrt(1+1.25*t))+140*15*int((z-k3)*normpdf(z,0,1),z,k3,inf)/(m3*t);
t0=sqrt((100+110+120)/(60*11.5+50*14+40*9));
tmin=1e-2;
tmax=4/3;
k1min=0;
k2min=0;
k3min=0;
k1max=2;
k2max=2;
k3max=2;

g1=1-200*sqrt(1+1.25*t)*int((z-k1)*normpdf(z,0,1),z,k1,inf)/(50*t);%%将约束标准化，因为3个约束常数相同，故不必标准化
g2=1-240*sqrt(1+1.25*t)*int((z-k2)*normpdf(z,0,1),z,k2,inf)/(50*t);  
g3=1-300*sqrt(1+1.25*t)*int((z-k3)*normpdf(z,0,1),z,k3,inf)/(50*t);

NP=50; %进化多少代

%%初始化种群，种群长度40
size=40;
E=zeros(40,7); %前四列为初始解，第五列为适应函数值，第六列记录是否为可行解,第七列记录约束函数差值
E(:,1)=tmin+(tmax-tmin)*rand(size,1);
E(:,2)=k1min+(k1max-k1min)*rand(size,1);
E(:,3)=k2min+(k2max-k2min)*rand(size,1);
E(:,4)=k3min+(k3max-k3min)*rand(size,1);
fv=inf;%初始最优值为无穷大的值
D=zeros(NP,4);%用来记录每代的最优解,平均值，最差解,最优解是否为可行解

%%计算适应函数罚函数值，判断是否为可行解
for i=1:size
    B=zeros(1,6);
    B(1)=subs(g1,[t,k1,k2,k3],E(i,(1:4)));
    B(2)=subs(g2,[t,k1,k2,k3],E(i,(1:4)));
    B(3)=subs(g3,[t,k1,k2,k3],E(i,(1:4)));
    if min(B)>=0
       E(i,6)=1;
       E(i,5)=subs(f,[t,k1,k2,k3],E(i,(1:4)));
    else
       E(i,6)=0;
       E(i,5)=0;
    end
for j=1:3
    if B(j)>=0
       B(j)=0;
    else
       B(j)=abs(B(j));
    end
end
E(i,7)=sum(B);
end
fmax=max(E(:,5));
for i=1:size
if E(i,6)<1e-6
   E(i,5)=fmax+E(i,7);
   end
end
  
%%遗传进化
for g=1:NP 
    %%竞标赛选择  %%小生态技术
    M=zeros(size,4);%用来存储优胜者的中间矩阵
        for i=1:size
            %A=randperm(size,6);
            %dij1=sqrt(0.5*(E(A(1),1)-E(A(2),1))^2);  %%小生态技术
            %dij2=sqrt(0.5*(E(A(1),1)-E(A(3),1))^2); 
            %dij3=sqrt(0.5*(E(A(1),1)-E(A(4),1))^2);
            %dij4=sqrt(0.5*(E(A(1),1)-E(A(5),1))^2);
            %dij5=sqrt(0.5*(E(A(1),1)-E(A(6),1))^2);
            %if dij1<0.1
                %if E(A(1),5)<=E(A(2),5)
                    %M(i,:)=E(A(1),(1:4));
                %else
                    %M(i,:)=E(A(2),(1:4));   
                %end
                %continue;
            %elseif dij2<0.1
                %if E(A(1),5)<=E(A(3),5)
                    %M(i,:)=E(A(1),(1:4));
                %else
                    %M(i,:)=E(A(3),(1:4));   
                %end
                %continue;
            %elseif dij3<0.1
                %if E(A(1),5)<=E(A(4),5)
                    %M(i,:)=E(A(1),(1:4));
                %else
                    %M(i,:)=E(A(4),(1:4));   
                %end
                %continue;
            %elseif dij4<0.1
                %if E(A(1),5)<=E(A(5),5)
                    %M(i,:)=E(A(1),(1:4));
                %else
                    %M(i,:)=E(A(5),(1:4));   
                %end
                %continue;
            %elseif dij5<0.1
                %if E(A(1),5)<=E(A(6),5)
                    %M(i,:)=E(A(1),(1:4));
                %else
                    %M(i,:)=E(A(6),(1:4));   
                %end
            %else
                %M(i,:)=E(A(1),(1:4));
            %end  
            A=randperm(size,2);
            if E(A(1),5)<=E(A(2),5)
                M(i,:)=E(A(1),(1:4));
            else
                M(i,:)=E(A(2),(1:4));
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
    %%只在可行解时出错是怎么回事？是不是变异的原因
            if rand()>0.5
                c=rand();
                x2=max(M(A(1),2),M(A(2),2));
                x1=min(M(A(1),2),M(A(2),2));
                beita1_t=1+2*(x1-k1min)/(x2-x1);
                rfa_t=2-beita1_t^(-2);
                if c<=1/rfa_t
                    beita2_t=sqrt(rfa_t*c);
                else
                    beita2_t=sqrt(1/(2-rfa_t*c));
                end
                E(2*j-1,2)=0.5*(x1+x2-beita2_t*(x2-x1));
                E(2*j,2)=0.5*(x1+x2+beita2_t*(x2-x1));
            end
            
            if rand()>=0.5
                c=rand();
                x2=max(M(A(1),3),M(A(2),3));
                x1=min(M(A(1),3),M(A(2),3));
                beita1_t=1+2*(x1-k2min)/(x2-x1);
                rfa_t=2-beita1_t^(-2);
                if c<=1/rfa_t
                    beita2_t=sqrt(rfa_t*c);
                else
                    beita2_t=sqrt(1/(2-rfa_t*c));
                end
                E(2*j-1,3)=0.5*(x1+x2-beita2_t*(x2-x1));
                E(2*j,3)=0.5*(x1+x2+beita2_t*(x2-x1));
            end
            
            if rand()>0.5
                c=rand();
                x2=max(M(A(1),2),M(A(2),4));
                x1=min(M(A(1),2),M(A(2),4));
                beita1_t=1+2*(x1-k3min)/(x2-x1);
                rfa_t=2-beita1_t^(-2);
                if c<=1/rfa_t
                    beita2_t=sqrt(rfa_t*c);
                else
                    beita2_t=sqrt(1/(2-rfa_t*c));
                end
                E(2*j-1,4)=0.5*(x1+x2-beita2_t*(x2-x1));
                E(2*j,4)=0.5*(x1+x2+beita2_t*(x2-x1));
            end
        end

        %%变异
         for i=1:size
            nita=100+g;
            pm=1/size+g*(1-1/size)/NP;
            if rand()<pm
                u=rand();
                x=E(i,1);
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
                   delta_2=(2*u)^(1/(nita+1))-1;
                else
                   delta_2=1-(2*(1-u))^(1/(nita+1));
                end
                if x+delta_2*deltamax>=k1min
                    E(i,2)=x+delta_2*deltamax;
                end
            end
            
            if rand()<pm
                u=rand();
                x=E(i,3);
                deltamax=1;
                if u<=0.5
                   delta_2=(2*u)^(1/(nita+1))-1;
                else
                   delta_2=1-(2*(1-u))^(1/(nita+1));
                end
                if x+delta_2*deltamax>=k2min
                    E(i,3)=x+delta_2*deltamax;
                end
            end
            
             if rand()<pm
                u=rand();
                x=E(i,4);
                deltamax=1;
                if u<=0.5
                   delta_2=(2*u)^(1/(nita+1))-1;
                else
                   delta_2=1-(2*(1-u))^(1/(nita+1));
                end
                if x+delta_2*deltamax>=k3min
                    E(i,4)=x+delta_2*deltamax;
                end
             end  
         end
             
        for i=1:size
    %%计算子代的适应罚函数值，判断是否为可行解
            B=zeros(1,6);
            B(1)=subs(g1,[t,k1,k2,k3],E(i,(1:4)));
            B(2)=subs(g2,[t,k1,k2,k3],E(i,(1:4)));
            B(3)=subs(g3,[t,k1,k2,k3],E(i,(1:4)));
            if min(B)>=0
                E(i,6)=1;
                E(i,5)=subs(f,[t,k1,k2,k3],E(i,(1:4)));
            else
                E(i,6)=0;
                E(i,5)=0;
            end
            for j=1:3
                if B(j)>=0
                    B(j)=0;
                else
                    B(j)=abs(B(j));
                end
            end
            E(i,7)=sum(B);
        end
    fmax=max(E(:,5));
    for i=1:size
        if E(i,6)<1e-6
            E(i,5)=fmax+E(i,7);
        end
    end
    [Q,IX]=sort(E,1);
    D(g,1)=Q(1,5);
    D(g,2)=mean(Q(:,5));
    D(g,3)=Q(size,5); 
    D(g,4)=E(IX(1,5),6);
    if Q(1,5)<fv && D(g,4)==1
        fv=Q(1,5);
        xm=E(IX(1,5),(1:4));
        xm=[xm,E(IX(1,5),6)];
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
legend('最优值','平均值','最差值');  %%画图时只显示出可行解
hold off
D
end

 