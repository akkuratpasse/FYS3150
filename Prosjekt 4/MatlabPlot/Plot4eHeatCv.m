% %Disorder
% t = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/t_Spins20.txt');
% HeatCv20 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/HeatCv20.txt');
% HeatCv40 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/HeatCv40.txt');
% HeatCv60 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/HeatCv60.txt');
% HeatCv80 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/HeatCv80.txt');
% 
% 
% plot(t, HeatCv20, 'r')
% hold on
% plot(t, HeatCv40, 'g')
% hold on 
% plot(t, HeatCv60, 'b')
% hold on
% plot(t, HeatCv80, 'y')
% hold off
% title('HeatCv: mcs = 10000, n temp = 20, Disorder')
% xlabel('k_bT') % x-axis label
% ylabel('HeatCv') % y-axis label
% legend('HeatCv20', 'HeatCv40', 'HeatCv60', 'HeatCv80')

%Order
t = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/t.txt');
HeatCv20 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/HeatCv20.txt');
HeatCv40 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/HeatCv40.txt');
HeatCv60 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/HeatCv60.txt');
HeatCv80 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/HeatCv80.txt');


plot(t, HeatCv20, 'r')
hold on
plot(t, HeatCv40, 'g')
hold on 
plot(t, HeatCv60, 'b')
hold on
plot(t, HeatCv80, 'y')
hold off
title('HeatCv: mcs = 10000, n temp = 20, ORDER')
xlabel('k_bT') % x-axis label
ylabel('HeatCv') % y-axis label
legend('HeatCv20', 'HeatCv40', 'HeatCv60', 'HeatCv80')