function Draw3DBar

x = 0:10;
B = 23:2:49;
Q = [0	0	7	9	11	13	15	17	18	18	18	18	18	18
0	0	7	9	11	13	15	17	17	17	17	17	17	17
0	0	7	9	11	13	15	16	16	16	16	16	16	16
0	0	7	9	11	13	15	15	15	15	15	15	15	15
0	0	0	9	11	13	14	14	14	14	14	14	14	14
0	0	0	9	11	13	13	13	13	13	13	13	13	13
0	0	0	9	11	12	12	12	12	12	12	12	12	12
0	0	0	9	11	11	11	11	11	11	11	11	11	11
0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0
];

bar3(Q);
set(gca,'XTickLabel',B)
set(gca,'YTickLabel',x)
xlabel('R', 'FontWeight', 'bold');
ylabel('x', 'FontWeight', 'bold');
zlabel('Q^\ast', 'FontWeight', 'bold');

% figure()
% y = Q;
% for i = 1 : size(Q,1)
%     y(i, :) = Q(i, :) + x(i);
% end
% bar3(y);
% set(gca,'XTickLabel',B)
% set(gca,'YTickLabel',x)
% xlabel('B');
% ylabel('x');
% zlabel('y^\ast');

end