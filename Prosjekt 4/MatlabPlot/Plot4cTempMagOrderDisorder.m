x = linspace(0,10000,10000);
% ORDER
OMT1 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/OrdermagnetizationsT=1.000000.txt');
OMT24 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/OrdermagnetizationsT=2.400000.txt');
% DISORDER
DMT1 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/DISORDERmagnetizationsT=1.000000.txt');
DMT24 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/DISORDERmagnetizationsT=2.400000.txt');

%Zoom
figure
subplot(2,2,1);
plot(x,OMT1/(20))
hold on
plot(x,DMT1/(20))
title('Temperature = 1, nspins = 20')
xlabel('#Mc') % x-axis label
ylabel('\langle M\rangle') % y-axis label
legend('OMT1 = Order','DMT1 = Disorder ')
hold off

subplot(2,2,2);
plot(x,OMT1/(20))
hold on
plot(x,DMT1/(20))
title('Temperature = 1, nspins = 20')
xlabel('#Mc') % x-axis label
ylabel('\langle M\rangle') % y-axis label
legend('OMT1 = Order','DMT1 = Disorder ')
hold off

subplot(2,2,3);
plot(x,OMT24/(20))
hold on
plot(x,DMT24/(20))
title('Temperature = 2.4, nspins = 20')
xlabel('#Mc') % x-axis label
ylabel('\langle M\rangle') % y-axis label
legend('OMT24 = Order','DMT24 = Disorder ')
hold off

subplot(2,2,4);
plot(x,OMT24/(20))
hold on
plot(x,DMT24/(20))
title('Temperature = 2.4, nspins = 20')
xlabel('#Mc') % x-axis label
ylabel('\langle M\rangle') % y-axis label
legend('OMT24 = Order','DMT24 = Disorder ')
hold off