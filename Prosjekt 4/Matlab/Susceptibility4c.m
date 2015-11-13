%Temperature
t = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4c-Desktop_Qt_5_5_0_clang_64bit-Debug/t.txt');
%Energy 
AE = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4c-Desktop_Qt_5_5_0_clang_64bit-Debug/AE.txt');
%Magnetization
AabsM = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4c-Desktop_Qt_5_5_0_clang_64bit-Debug/AabsM.txt');

figure
subplot(2,1,1);
plot(t, AE, 'r')
title('Susceptibility: mcs = 10000, nspins = 40, ntemp = 20')
xlabel('Temperature') % x-axis label
ylabel('AE') % y-axis label
subplot(2,1,2);
plot(t, AabsM, 'b')
title('HeatCv: mcs = 10000, nspins = 40, ntemp = 20')
xlabel('Temperature') % x-axis label
ylabel('AabsM') % y-axis label
