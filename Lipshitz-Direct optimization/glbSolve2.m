function [x_k1,glbfMin,safeinventory,reorder,leadtime,t,iter] = glbSolve2(x_L,x_U,EpsGlob,MaxFunc,var)
format long;
% clear;
if nargin~=5
   disp('the wrong input');
end

% syms k q;
%[x_k1,glbfMin,safeinentory,reorder,leadtime,t,iter] = glbSolve2([0;0],[100;10000],1e-6,100,[k,q])
%MaxFunc为最大迭代次数

DD=25;sigma=10;A=55;h=4;lc=230;aa=9.8;bb=0.001;cc=0.03;
syms k q z;
intloss=matlabFunction(int((z-k)*normpdf(z,0,1),z,k,inf)); %定义缺货积分函数
f=A*DD/q+h*(q/2+k*sigma*sqrt(aa/q+bb*q+cc))+lc*DD*sigma*sqrt(aa/q+bb*q+cc)*intloss(k)/q; %定义目标函数


% %初始值
% q0=sqrt(2*A*D/h);
% k0=norminv(h*q0/lc/D);
% 
% %q上界
% q_upp=2*subs(f,var,[q0,k0])/h;
%  STEP 1, Initialization
   nFunc     = 1;
   x_D =x_U-x_L;  % x_L是两个变量的下届

   n   = length(x_L);  % Problem dimension问题的维度
 
   tol2=1e-8;
   m = 1;             % Current number of rectangles 目前的矩形数目
   mxDim = min(MaxFunc+n*209,999); %问题最大维度？
   C = [0.5*ones(n,1), zeros(n,mxDim)]; %中心点
   
   x_m = x_L + 0.5*x_D;   % Start with mid point  初始解
   glbfMin =subs(f,var,x_m);  % Function value at x_m 初始值
 
   iMin = 1; % The rectangle which minimizes (F - glbfMin + E)./D where
             % E = max(EpsGlob*abs(glbfMin),1E-8)

% Matrix with all rectangle side lengths in each dimension
   L = [0.5*ones(n,1), zeros(n,mxDim)];  %L用来记录矩形的中心和边长
   % Vector with distances from centerpoint to the vertices
   %D = sqrt(sum(L.^2));  
   D = [sqrt(n/4),zeros(1,mxDim)];  %D用来记录中心点到顶点的距离
   % Vector with function values
   F = [glbfMin,zeros(1,mxDim)];   %F用来记录目标函数值

   DSort = D(1);     % Row vector of all different distances, sorted
   DMin  = glbfMin(1);  % Row vector of minimum function value for each distance   
   c23 = 2/3;
   iter   = 0; 
%    Iter <= MaxIter  
t1=clock;

S = [];S_2 = [];  idx = find(DSort==D(iMin)); %最小值可能存在的搜索区域S
   S_1 = zeros(1,length(DSort)-idx+1);
% STEP2 ITERATION LOOP&&nFunc <= MaxFunc
flag1=1;
counter=0;
iter=1;
temp=10000;
while  iter<=26
  
   %  STEP 2  Identify the set S of all potentially optimal rectangles
     % Set of all potentially optimal rectangles
  
   mz = 0;
   for i = idx : length(DSort)
      idx2 = find(  F(1:m)==DMin(i)  &  D(1:m)==DSort(i)  );
      S_1(mz+1:mz+length(idx2)) = idx2;
      mz = mz + length(idx2);
   end
   
   if length(DSort)-idx > 1
      a1 = D(iMin);
      b1 = F(iMin);
      a2 = DSort(length(DSort));
      b2 = DMin(length(DSort));
      % The line is defined by: y = slope*x + const
      slope = (b2-b1)/(a2-a1);
      const = b1 - slope*a1;
      S_2 = iMin;
      for i = 1 : length(S_1)
         j = S_1(i);
         if j ~= iMin
            if F(j) <= slope*D(j) + const + tol2 
               S_2 = [S_2 j];
            end   
         end
      end
      % S_2 now contains all points in S_1 which lies on or below the line
        
      % Find the points on the convex hull defined by the points in S_2
      xx = D(S_2);
      yy = F(S_2);
      h  = convhull(xx,yy); % conhull is an internal subfunction
%   h  = tomsol(11,xx,yy); % conhull is an internal subfunction
      S_3 = S_2(h);
   else
      S_3 = S_1;
   end
   S = S_3; 
   
   
   %  STEP 3, 5  Select any rectangle j in S
   for jj = 1:length(S)  
      j = S(jj);
      [max_L, I] = max(L(:,j));
      I = find( L(:,j)==max_L );
      delta = c23*max_L;
      w=zeros(length(I),1);
      mz=m;
      if m+2*length(I) > size(C,2)
         C  = [C, zeros(n,MaxFunc)];
         F  = [F, zeros(1,MaxFunc)];
         L  = [L, zeros(n,MaxFunc)];
      end
      for ii = 1:length(I) 
         i = I(ii);
         e_i=[zeros(i-1,1);zeros(n-i,1)];
         c_m1=C(:,j)+delta*e_i;
%          c_m1 = C(:,j);
         c_m1(i)=c_m1(i)+delta;
         x_m1 = x_L+c_m1.*x_D; 
         f_m1 = subs(f,[var(1),var(2)],x_m1);
         c_m2 = C(:,j);  
         c_m2(i)=c_m2(i)-delta;
         x_m2 = x_L+ c_m2.*x_D;
         f_m2 = subs(f,[var(1),var(2)],x_m2);
         nFunc=nFunc+2; 
         if abs(f_m2)<=abs(f_m1)
         w(ii) = f_m2;
%          x_o=x_m2;
         else
             w(ii)=f_m1;
%          x_o=x_m1;
         end
%          x_o1=x_m1;x_o2=x_m2;
%          x1=x_o1(1);y1=x_o1(2);
%          x2=x_o2(1);y2=x_o2(2);
%          plot(x1,y1,'+r');AXIS([0 1.1 0 1.1]);set(gca, 'YGrid', 'on');set(gca, 'XGrid', 'on');
%          hold on;  
%          plot(x2,y2,'+r');AXIS([0 1.1 0 1.1]);set(gca, 'YGrid', 'on');set(gca, 'XGrid', 'on');
%          hold on;
         mz = mz+2;
         C(:,mz-1:mz) = [c_m1 c_m2];
         F(:,mz-1:mz) = [f_m1 f_m2];          
      end
      [a b] = sort(w);
      glbfMin = min(glbfMin,a(1)); 
      for ii = 1:length(I)
         i = I(b(ii));         
         ix1 = m + 2*b(ii)-1;       
         ix2 = ix1 + 1;         
         L(i,j) = delta/2;
         L(:,[ix1,ix2]) = L(:,[j,j]);
         D([j ix1 ix2])   = norm(L(:,j));
      end
      m = m + 2*length(I);
   end
   
   % UPDATE:
   glbfMin = min(F(1:m));  

   E = max(EpsGlob*abs(glbfMin),1E-8);
   DSort = D(1:m); % Use DSort instead of D before sorting in min computation
   [dummy iMin] = min( (F(1:m) - glbfMin + E)./DSort );
 
   i = 1;
   while 1
      d_tmp = DSort(i);
      idx = find(DSort~=d_tmp);
      DSort = [d_tmp DSort(idx)];
      if i==length(DSort)
         break;
      else
         i = i + 1;
      end
   end
   DSort = sort(DSort);   
   DMin = zeros(1,length(DSort));
   for i = 1:length(DSort);
      idx1 = find(D(1:m)==DSort(i));
      %idx1 = find( abs( D-DSort(i) ) <= tol1 );
      DMin(i) = min(F(idx1));
   end
   f_m5=glbfMin;
    if (abs(temp-f_m5)/temp)>=1e-6
        counter=0;
    else
        counter=counter+1;
        if counter>10
            flag1=0;
        end
    end
     iter = iter + 1; 
   disp('iter=');disp(iter);
%    disp('abs(temp-f_m5)=');disp(abs(temp-f_m5));
     temp=f_m5;
%      disp('temp');disp(temp);
end % ITERATION LOOP

C=C(:,1:m);F=F(:,1:m);
B=find(F==glbfMin);
    x_k1= x_L+C(:,B(1)).*x_D; 
t2=clock;
t=etime(t2,t1);
safeinventory=x_k1(1)*sigma*sqrt(aa/x_k1(2)+bb*x_k1(2)+cc);
leadtime=aa/x_k1(2)+bb*x_k1(2)+cc;
reorder=DD*(aa/x_k1(2)+bb*x_k1(2))+safeinventory;