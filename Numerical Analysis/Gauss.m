function x=Gauss(A,b)  
    n=size(A,1);
    x=zeros(n,1);
    for k=1:n-1
        if A(k,k)==0
            error('algorithm failed');
        end
        for i=(k+1):n
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