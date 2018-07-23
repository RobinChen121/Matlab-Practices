function NewtonOptim(f,x0,eps)
X=symvar(f);
df1=jacobian(f,X);
df2=jacobian(df1,X);
mark=eval(norm(subs(df1,X,x0),2));
while mark>eps
    d=-subs(df1,X,x0)/(subs(df2,X,x0));
    x0=x0+d;
    mark=eval(norm(subs(df1,X,x0),2));
end
fval=subs(f,X,x0);


end