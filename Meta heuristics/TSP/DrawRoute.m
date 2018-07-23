function m=DrawRoute(Clist,BSF,bsf,p)
CityNum=size(Clist,1);
for i=1:CityNum-1
    plot([Clist(BSF(i),1),Clist(BSF(i+1),1)],[Clist(BSF(i),2),Clist(BSF(i+1),2)],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g');%用直线将城市连起来
    hold on;  % 允许在同一坐标系绘制不同的图形
end
plot([Clist(BSF(CityNum),1),Clist(BSF(1),1)],[Clist(BSF(CityNum),2),Clist(BSF(1),2)],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g');
title([num2str(CityNum),'城市TSP']);
text(5,5,['第 ',int2str(p),' 步','  最短距离为 ',num2str(bsf)]);
hold off;
pause(0.05); %matlab执行命令间隔时间，能够显示动态效果
