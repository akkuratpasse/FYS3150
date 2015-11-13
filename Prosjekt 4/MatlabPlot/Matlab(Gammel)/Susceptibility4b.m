t = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV3-Desktop_Qt_5_5_0_clang_64bit-Debug/t.txt');  % load temp
Susceptibility = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV3-Desktop_Qt_5_5_0_clang_64bit-Debug/Susceptibility.txt');
HeatCv = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Debug/HeatCv.txt');


figure
subplot(2,1,1);
plot(t, Susceptibility, 'r')
title('Susceptibility: mcs = 10000, n spins = 40, n temp = 20')
xlabel('T') % x-axis label
ylabel('Susceptibility') % y-axis label
subplot(2,1,2);
plot(t, HeatCv, 'b')
title('HeatCv: mcs = 10000, n spins = 40, n temp = 20')
xlabel('T') % x-axis label
ylabel('HeatCv') % y-axis label
