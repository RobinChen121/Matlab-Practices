function drawTSP10(Clist,BSF,qs,bsf,p,f)
CityNum=size(Clist,1)-1;
CarNum=size(qs,2)-1;
%col='gbymcwkgbymcwkgbymcwk';%让不同车辆的行驶路线颜色不一样
%line_s=struct('f',{'-','-','-','--','-',':','--','-.'});
line_c=[0.8,2,3.5,0.8,2,3.5,0.8,2,3.5,0.8,2,3.5,0.8,2,3.5];
scatter(Clist(1:CityNum,1),Clist(1:CityNum,2),'k');
hold on;
for i=1:CityNum
    text(Clist(i,1)+2,Clist(i,2)+2,num2str(i));  %画图错了，计算没错
end
for j=1:CarNum
    %annotation('arrow',[Clist(CityNum+1,1),Clist(BSF(qs(j)),1)]/200,[Clist(CityNum+1,2),Clist(BSF(qs(j)),2)]/200,'LineStyle','-');
    %text([Clist(CityNum+1,1),Clist(BSF(qs(j)),1)],[Clist(CityNum+1,2),Clist(BSF(qs(j)),2)],'chen zhen');
    plot([Clist(CityNum+1,1),Clist(BSF(qs(j)),1)],[Clist(CityNum+1,2),Clist(BSF(qs(j)),2)],strcat('k',':'),'LineWidth',line_c(j))
    hold on;
    for i=qs(j):qs(j+1)-2
        %annotation('arrow',[Clist(BSF(i),1),Clist(BSF(i+1),1)]/200,[Clist(BSF(i),2),Clist(BSF(i+1),2)]/200,'LineStyle','-')
        plot([Clist(BSF(i),1),Clist(BSF(i+1),1)],[Clist(BSF(i),2),Clist(BSF(i+1),2)],strcat('k','-'),'LineWidth',line_c(j));
        hold on;
    end
    %annotation('arrow',[Clist(BSF(qs(j+1)-1),1),Clist(CityNum+1,1)]/200,[Clist(BSF(qs(j+1)-1),2),Clist(CityNum+1,2)]/200,'LineStyle','-')
    plot([Clist(BSF(qs(j+1)-1),1),Clist(CityNum+1,1)],[Clist(BSF(qs(j+1)-1),2),Clist(CityNum+1,2)],strcat('k','-'),'LineWidth',line_c(j));
end
hold on;
plot(Clist(CityNum+1,1),Clist(CityNum+1,2),'ks','MarkerSize',13);
%axis([0,1,0,1]);
title([num2str(CityNum),'城市TSP']);
%legend('第一辆车','第二辆车','第三辆车');
hold off;
%if f==0
   % text(1,1,['第 ',int2str(p),' 步','  最小费用为 ',num2str(bsf,5),' 元']);
%else
   % text(1,1,['最终搜索结果：最小费用 ',num2str(bsf,5),' 元']);
%end
pause(0.05); 