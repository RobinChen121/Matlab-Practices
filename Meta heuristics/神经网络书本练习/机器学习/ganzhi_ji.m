function ganzhi_ji
    [a1,a2,a3]=textread('dataset1.txt','%f%f%c'); %%读取数据
    [b1,b2,b3]=textread('dataset2.txt','%f%f%c');
    x1=[a1;b1];
    x2=[a2;b2];
    b=[a3;b3];
    
    %%将男女类别变为+-1
    n=size(x1,1);
    y=zeros(n,1);
    index1=0;
    index2=0;
    for i=1:n
        if b(i)=='m'||b(i)=='M'
            y(i)=1;
            index1=index1+1;
            p1(:,index1)=[x1(i);x2(i)];
            t1(index1)=1;
        elseif b(i)=='f'||b(i)=='F'
                y(i)=-1;
                index2=index2+1;
                p2(:,index2)=[x1(i);x2(i)];
                t2(index2)=0;
        end
    end
    p=[p1 p2];
    t=[t1 t2];
    plotpv(p,t);
    net=newp([140 200;20 140],1);
    [net,v,e]=adapt(net,p,t);
    [c1,c2,c3]=textread('dataset1.txt','%f%f%c');
    c=[c1 c2]';
    tt=sim(net,c);
    
end
    
    
    
    %for i=1:index1
        %plot(p1(1,i),p1(2,i),'o');
   % end
    %hold on
    %for i=1:index2
        %plot(p2(1,i),p2(2,i),'+');
    %end
    
    
    %%初始化
    %w0_1=1;
    %w0_2=1;
    %b0=1;
    %nita=0.8;
    
    %for i=1:n
        %if y(i)*(w0_1*x1(i)+w0_2*x2(i)+b0)<=0
            %b0=b0+nita*y(i);
            %w0_1=w0_1+nita*y(i)*x1(i);
            %w0_2=w0_2+nita*y(i)*x2(i);
        %end
    %end
    
    %%测试第三组数据
    %[c1,c2,c3]=textread('dataset1.txt','%n%n%c');
    %n2=size(c1,1);
    %y_t=zeros(n2,1);
    %y_c=zeros(n2,1);
    %for i=1:n2
        %if c3(i)=='m'||c3(i)=='M'
            %y_t(i)=1;
        %end
        %elseif c3(i)=='f'||c3(i)=='F'
                %y_t(i)=-1;
        %end
    %end
    
    %%验证测试数据
    %index=0;
    %for i=1:n2
        %z=w0_1*c1(i)+w0_2*c2(i)+b0;
        %if z>0
            %y_c(i)=1;
        %else
            %y_c(i)=-1;
        %end
        %if y_c(i)==y_t(i)
            %index=index+1;
        %end
    %end
    %disp(index/n2);
%end