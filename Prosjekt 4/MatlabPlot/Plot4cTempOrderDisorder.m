
x = linspace(0,10000,10000);
% ORDER
OT1 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/ORDERenergiesT=1.000000.txt');
OT24 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/ORDERenergiesT=2.400000.txt');
% DISORDER
DT1 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/DISORDERenergiesT=1.000000.txt');
DT24 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/DISORDERenergiesT=2.400000.txt');


% figure
% subplot(2,1,1);
% plot(x,OT1/(40^2))
% hold on
% plot(x,DT1/(40^2))
% title('Temperature = 1, nspins = 20')
% xlabel('#Mc') % x-axis label
% ylabel('\langle E\rangle') % y-axis label
% legend('OT1 = Order','DT1 = Disorder ')
% hold off
% 
% subplot(2,1,2);
% plot(x,OT24/(40^2))
% hold on
% plot(x,DT24/(40^2))
% title('Temperature = 2.4, nspins = 20')
% xlabel('#Mc') % x-axis label
% ylabel('\langle E\rangle') % y-axis label
% legend('OT24 = Order','DT24 = Disorder ')
% hold off

%Zoom
figure
subplot(2,2,1);
plot(x,OT1/(40^2))
hold on
plot(x,DT1/(40^2))
title('Temperature = 1, nspins = 20')
xlabel('#Mc') % x-axis label
ylabel('\langle E\rangle') % y-axis label
legend('OT1 = Order','DT1 = Disorder ')
hold off

subplot(2,2,2);
plot(x,OT1/(40^2))
hold on
plot(x,DT1/(40^2))
title('Temperature = 1, nspins = 20')
xlabel('#Mc') % x-axis label
ylabel('\langle E\rangle') % y-axis label
legend('OT1 = Order','DT1 = Disorder ')
hold off

subplot(2,2,3);
plot(x,OT24/(40^2))
hold on
plot(x,DT24/(40^2))
title('Temperature = 2.4, nspins = 20')
xlabel('#Mc') % x-axis label
ylabel('\langle E\rangle') % y-axis label
legend('OT24 = Order','DT24 = Disorder ')
hold off

subplot(2,2,4);
plot(x,OT24/(40^2))
hold on
plot(x,DT24/(40^2))
title('Temperature = 2.4, nspins = 20')
xlabel('#Mc') % x-axis label
ylabel('\langle E\rangle') % y-axis label
legend('OT24 = Order','DT24 = Disorder ')
hold off