function x=lie_Gauss(A,b)
    n=size(A,1);
    m=size(A,2);
    for k=1:n-1
        max=abs(A(k,k));
        i_max=k;
        for i=k+1:n
            if abs(A(i,k))>abs(A(i-1,k))
                max=abs(A(i,k));
                i_max=i;
            end
        end
        for j=1:m
            temp=A(k,j);
            A(k,j)=A(i_max,j);
            A(i_max,j)=temp;
            temp=b(k);
            b(k)=b(i_max);
            b(i_max)=temp;
        end
        for i=k+1:n
            m=A(i,k)/A(k,k);
            for j=(k+1):n
                A(i,j)=A(i,j)-m*A(k,j);  
            end
            for j=1:k
                A(i,j)=0;  
            end
            b(i)=b(i)-m*b(k);
        end
    end
     x(n)=b(n)/A(n,n);
    for k=n-1:-1:1
        temp=0.0;
        for j=k+1:n
            temp=temp+A(k,j)*x(j);
        end
        x(k)=(b(k)-temp)/A(k,k);
    end
end