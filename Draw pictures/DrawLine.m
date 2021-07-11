function  DrawLine

x = [10000, 15000, 20000, 25000, 30000];
y = [6.51, 0.34, 0.14, 0.09, 0.04];
p = plot(x,y, 'ro-');
a = ancestor(p, 'axes');
a.YAxis.TickLabelFormat = '%g%%';
ylim([-1, 8]);
set(gca,'xtick',0:5000:35000)
ylabel('Profit gap');
xlabel('Initial cash $C_0$','Interpreter','latex');

end

