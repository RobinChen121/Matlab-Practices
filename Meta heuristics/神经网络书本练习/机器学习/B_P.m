%clc;
%clear;
%载入数据
[a1,a2,a3]=textread('dataset1.txt','%f%f%c');
[b1,b2,b3]=textread('dataset1.txt','%f%f%c');
[c1,c2,c3]=textread('dataset1.txt','%f%f%c');

%将性别用数字矩阵表示出来
n1=length(a3);
n2=length(b3);
n3=length(c3);
output=zeros(n1+n2+n3,2);
for i=1:n1
    switch a3(i)
        case 'm'
            output(i,:)=[1 0];
        case 'M'
            output(i,:)=[1 0];
        case 'f'
            output(i,:)=[0 1];
        case 'F'
            output(i,:)=[0 1];
    end
end
for i=1:n2
    switch b3(i)
        case 'm'
            output(n1+i,:)=[1 0];
        case 'M'
            output(n1+i,:)=[1 0];
        case 'f'
            output(n1+i,:)=[0 1];
        case 'F'
            output(i,:)=[0 1];
    end
end
for i=1:n3
    switch c3(i)
        case 'm'
            output(n1+n2+i,:)=[1 0];
        case 'M'
            output(n1+n2+i,:)=[1 0];
        case 'f'
            output(n1+n2+i,:)=[0 1];
        case 'F'
            output(i,:)=[0 1];
    end
end

%前两组数据作为训练数据，最后一组数据作为测试数据
input_train=[a1 a2;b1 b2]';
input_test=[c1 c2]';
output_train=output(1:n1+n2,:)';
output_test=output(n1+n2+1:n1+n2+n3,:)';

%输入数据归一化
inputn= input_train./repmat(sqrt(sum(input_train.*input_train)),size(input_train,1),1);

%网络结构
innum=2;
midnum=3;
outnum=2;
nita=0.8;

%权值阀值初始化
w1=rands(midnum,innum);
b1=rands(midnum,1);
w2=rands(midnum,outnum);
b2=rands(outnum,1);

E=zeros(20,1); %计算每次的总误差
for ii=1:20  %训练20次
    E(ii)=0;
    for i=1:n1+n2
        x=inputn(:,i);  %选择本次训练数据
        
        %隐含层输出
        I=zeros(midnum,1);
        Iout=zeros(1,outnum);
        for j=1:midnum
            I(j)=x'*w1(j,:)'+b1(j);
            Iout(j)=1/(1+exp(-I(j)));
        end
        %输出层输出
        yn=w2'*Iout'+b2;
        
        %预测误差
        e=output_train(:,i)-yn;
        E(ii)=E(ii)+sum(abs(e));
        
        %计算w2,b2调整量
        dw2=e*Iout;
        db2=e';
        
        %计算w1,b1调整量
        FI=zeros(1,midnum);
        for j=1:midnum
            S=1/(1+exp(-I(j)));
            FI(j)=S*(1-S);
        end
        dw1=zeros(innum,midnum);
        db1=zeros(1,midnum);
        for k=1:innum
            for j=1:midnum
               dw1(k,j)=FI(j)*x(k)*(e(1)*w2(j,1)+e(2)*w2(j,2));    
               db1(j)=FI(j)*(e(1)*w2(j,1)+e(2)*w2(j,2));
            end
        end
        w1=w1+nita*dw1';
        b1=b1+nita*db1';
        w2=w2+nita*dw2';
        b2=b2+nita*db2';
    end    
end

%输入数据归一化
inputn_test= input_test./repmat(sqrt(sum(input_test.*input_test)),size(input_test,1),1);
fore=zeros(2,n3);
for i=1:n3
    for j=1:midnum
        I(j)=inputn_test(:,i)'*w1(j,:)'+b1(j);
        Iout(j)=1/(1+exp(-I(j)));
    end
    fore(:,i)=w2'*Iout'+b2;
end

%output_fore=zeros(n3,1);
for i=1:n3
    output_fore(i)=find(fore(:,i)==max(fore(:,i)));
end
error=output_fore-output_test';
disp;




