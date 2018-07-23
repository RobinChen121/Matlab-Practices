function h=conhull(x,y) %x,y为横纵坐标
if nargin~=2
    disp('the wrong input');
end
[newx,index] = sort(x);
newy = y(index);
x=newx;
y=newy;    %%按x重新排序
m=length(x);
h=ones(1,m);
for i=1:m
    h(i)=h(i)+i-1;
end
if m>3
    start=1;
    v=start;
    w=m;
    flag=0;
    while next(v,m)~=start||flag==0
        if next(v,m)==w
            flag=1;
        end
        a=v;b=next(v,m);c=next(next(v,m),m);
        A=[x(a) y(a) 1;x(b) y(b) 1;x(c) y(c) 1];
        if det(A)>=0
            leftturn=0;
        else
            leftturn=1;
        end
        if leftturn
            v=next(v,m);
        else
            j=next(v,m);
            x(j)=[];y(j)=[];h(j)=[];
            m=m-1;w=w-1;
            if v==1
                v=m;
            else
                v=v-1;
            end
        end
    end
end
end

function v=next(v,m)
if v<m
    v=v+1;
elseif v==m
    v=1;
end
end
