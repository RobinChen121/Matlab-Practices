function Knapsack
c=[2 2 6 5 4];v=[6 3 5 4 6];
f=zeros(5,11);x=zeros(5,11);xx=zeros(5,1);
for i=5:-1:1
    for S=0:10
        if i==5
            f(i,S+1)=v(i)*floor(S/c(i));
            x(i,S+1)=floor(S/c(i));
        else
            xMax=floor(S/c(i));
            ff=zeros(xMax+1,1);
            for k=0:xMax
                ff(k+1)=v(i)*k+f(i+1,S-c(i)*k+1);
            end
            [f(i,S+1),index]=max(ff);
            x(i,S+1)=index-1;
        end
    end
end
[optValue,index]=max(f(1,:));
xx(1)=x(1,index);
tempS=index;
fprintf('optimal solution:%d\n',optValue);
for i=2:5
    xx(i)=x(i,tempS-c(i-1)*xx(i-1));
    tempS=tempS-c(i-1)*xx(i-1);
end
for i=1:5
    fprintf('put %d item%d in the bag\n',xx(i),i);
end

end