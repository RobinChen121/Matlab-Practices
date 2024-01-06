function MyWW
% x表示各阶段生产决策，为0-1变量
% y表示各阶段生产的数量
% I为各阶段库存水平
% I0 initial inventory level

I0=0;
n=6;

% d=randi([5 25],1,n); 
% S=randi([20 100],1,n); 
% h=randi([5 20],1,n); 
% c=randi([10 25],1,n); 

d=[4, 10, 5, 2, 10, 5]; 
S= 75 * ones(1, n);
h=5*ones(1,n);
c=0*ones(1,n);

x=zeros(1,n);
y=zeros(1,n);
I=zeros(1,n);

C=1e4*ones(n,n); %成本矩阵
opt_cost=zeros(n,1);
initial_I=zeros(n,1);
for i=1:n
    if i>1
        opt_cost=min(C(:,i-1));
    end
    for j=i:n
        if I0>=sum(d(1:j))
            if i==1
                C(i,j)=h(1:j)*(I0-cumsum(d(i:j)))';
            else
                C(i,j)=opt_cost+h(i:j)*(I0-cumsum(d(i:j)))';
            end
        else
            initial_I(i)=I0;
            if i>1
                if I0>sum(d(1:i-1))
                    initial_I(i)=I0-sum(d(1:i-1));
                else
                    initial_I(i)=0;
                end
            end
            prod_amount=sum(d(i:j))-initial_I(i);
            temp_I=initial_I(i)+prod_amount;
            h_sum=h(i:j)*(temp_I-cumsum(d(i:j)))';
            if i>1
                C(i,j)=opt_cost+S(i)+h_sum+c(i)*prod_amount;
            else
                C(i,j)=S(i)+h_sum+c(i)*prod_amount;
            end
        end
   end
end

%back track
j=n;
while j>=1
    [~,index]=min(C(:,j));
    if I0<sum(d(1:index))
        x(index)=1;
        if i>1
            y(index)=sum(d(index:j))-initial_I(index);            
        else
            y(index)=sum(d(index:j))-initial_I(index);
        end
    end
    j=index-1;
end
fprintf('各阶段库存量为:\n');
for i=1:n
    if i==1
        I(i)=I0+y(i)-d(i);
    else
        I(i)=I(i-1)+y(i)-d(i);
    end
    fprintf('  %d',I(i));
end
fprintf('\n');
fprintf('各阶段生产量为:\n');
for i=1:n
    fprintf('  %d',y(i));
end
fprintf('\n');
fprintf('各阶段是否生产:\n');
for i=1:n
    fprintf('  %d',x(i));
end
fprintf('\n');

end

