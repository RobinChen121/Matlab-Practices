function GM11
% 灰色预测 GM(1, 1)模型，比较适合指数趋势的数据
%

close all;
Y0 = [8353,5091,3052]; % 原始数据
Y1 = cumsum(Y0);
X = Y0(2 : end)';
N = length(Y0);
B = zeros(N - 1, 2);
for i = 1 : N - 1
    B(i, 1) = -0.5 * (Y1(i) + Y1(i + 1));
    B(i, 2) = 1;
end
alB = inv(B'*B)*B'*X;
alpha = alB(1);
mu= alB(2);

Y2 = zeros(1, N + 1);
Y2(1) = Y0(1);
for i = 2 : N + 1
    Y2(i) = (Y0(1) - mu/alpha) * exp(-alpha * (i-1)) + mu/alpha;
end
foreValue = zeros(1, N + 1);
foreValue(1) = Y0(1);
for i = 2 : N + 1
    foreValue(i) = Y2(i) - Y2(i - 1);
end
fprintf('forecast value is %.2f\n', foreValue(end));

%绘制曲线图
t1 = 1 : N;
t2 = 1 : N + 1;

Y0= [8353,5091,3052]
foreValue = [8353, 8353, 7048, 5449];
plot(t1, Y0,'ro'); hold on;
plot(t2, foreValue, 'g-');
legend('实际数据','预测数据');
xlim([1 5]);
ylim([0 9000]);
grid on;
end