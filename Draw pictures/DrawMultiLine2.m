function DrawMultiLine2
close all;

C0 = [20000	30000	40000	50000	60000];
VSS	= [0	31.2	7265.702171	1871.642851	1032.433361];
%EVPI =	[0	5891.65069	7642.595927	3235.77252	3266.427018];

VSS2	= [5.585082782	564.1721003	1878.962628	5956.541951	2176.683484];
%EVPI =	[10940.24966	11848.32593	14156.29598	19085.95082	14745.17495];

plot(C0, VSS, '-o', C0, VSS2, '-*');
ylim([0, 7500]);
xlim([10000, 60000])
xlabel('��ʼ�ʽ��� C_0');
ylabel('VSS');
legend('��������', '���ô���');


figure;
L = [0	1	2	3	4];
VSS	= [1402.526203	407.9113905	31	357	1270.188999];
%EVPI =	[4712.97156	2510.581451	5891.65069	4595.726477	4849.130759];

VSS2	= [1195.784427	509.4086787	564.1721003	491.2492881	951.2930959];
%EVPI =	[15132.7505	14124.41285	11848.32593	12926.23069	15668.32426];

plot(L, VSS, '-o', L, VSS2, '-*');
ylim([0, 1500]);
xlabel('�տ��ӳ��� L');
ylabel('VSS');
legend('��������', '���ô���');

figure;
H = [0	1000	2000	3000	4000];
VSS	= [272.2427044	43.13409852	31.2	0	0];
%EVPI =	[349.6579796	97.31655857	5891.65069	15.94548653	20.34638305];
VSS2	= [1071.332518	816.9264365	564.1721003	324.6863268	92.22263261];
%EVPI =	[12698.57437	12331.69068	11848.32593	11378.26124	11061.13849];

plot(H, VSS, '-o', H, VSS2, '-*');
ylim([0, 1200]);
xlabel('��ӷ��� H');
ylabel('VSS');
legend('��������', '���ô���');

figure;
B = [1000	5000	10000	20000	30000];
VSS	= [6.758472973	15.82081978	31.2 44.09635176	41.81541775];
VSS2 = [6.756	170	564.1721003	560.9690669	917.0017671];

plot(B, VSS, '-o', B, VSS2, '-*');
ylim([0, 1200]);
xlabel('������');
ylabel('VSS');
legend('��������', '���ô���');

figure;
r = [0.1	1.5	5	7.5	10];
VSS	= [89.85688585	31.22658783	1.810444192	1.810521976	153.7820997];
VSS2 = [564.1721003	514.1721003	464.1721003	484.1721003	414.1721003];

plot(r, VSS, '-o', r, VSS2, '-*');
ylim([0, 1200]);
xlabel('��������');
ylabel('VSS');
legend('��������', '���ô���');

figure;
r = [1	2	3	4	5];
VSS	= [1.062618638	31.22658783	0	1	2];
VSS2 = [512.1482787	31.22658783	0	31.09119852	2.005959108];

plot(r, VSS, '-o', r, VSS2, '-*');
ylim([0, 800]);
xlabel('���󲨶�ģʽ');
ylabel('VSS');
legend('��������', '���ô���');
end