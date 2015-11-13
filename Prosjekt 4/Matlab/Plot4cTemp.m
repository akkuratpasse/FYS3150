% Low Temp

%T1 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4cTemp-Desktop_Qt_5_5_0_clang_64bit-Release/energiesT=1.000000.txt');
%T24 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4cTemp-Desktop_Qt_5_5_0_clang_64bit-Release/energiesT=2.400000.txt')
x = linspace(0,10000,10000);

% ORDER
OT1 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/ORDERenergiesT=1.000000.txt');
OT24 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/ORDERenergiesT=2.400000.txt');
% DISORDER
DT1 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/DISORDERenergiesT=1.000000.txt');
DT24 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/DISORDERenergiesT=2.400000.txt');


figure
subplot(2,1,1);
plot(x,OT1)
hold on
plot(x,DT1)
title('Temperature = 1')
xlabel('#Mc') % x-axis label
ylabel('\langle E\rangle') % y-axis label
legend('OT1 = Order','DT1 = Disorder ')
hold off

subplot(2,1,2);
plot(x,OT24)
hold on
plot(x,DT24)
title('Temperature = 2.4')
xlabel('#Mc') % x-axis label
ylabel('\langle E\rangle') % y-axis label
legend('OT24 = Order','DT24 = Disorder ')
hold off