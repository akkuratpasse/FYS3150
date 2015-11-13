%Disorder
% t = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/t_Spins20.txt');
% AE20 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/AE20.txt');
% AE40 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/AE40.txt');
% AE60 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/AE60.txt');
% AE80 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/AE80.txt');
% 
% 
% plot(t, AE20, 'r')
% hold on
% plot(t, AE40, 'g')
% hold on
% plot(t, AE60, 'b')
% hold on
% plot(t, AE80, 'y')
% hold off
% title('\langle E\rangle: mcs = 10000, ntemp = 20, Disorder')
% xlabel('K_bT') % x-axis label
% ylabel('\langle E\rangle') % y-axis label
% legend('AE20','AE40','AE60', 'AE80' )
% 
%Order
t = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/t.txt');
AE20 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/AE20.txt');
AE40 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/AE40.txt');
AE60 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/AE60.txt');
AE80 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/AE80.txt');

plot(t, AE20/(20^2), 'r')
hold on
plot(t, AE40/(40^2), 'g')
hold on
plot(t, AE60/(60^2), 'b')
hold on
plot(t, AE80/(80^2), 'y')
hold off
title('\langle E\rangle: mcs = 10000, ntemp = 20, ORDER')
xlabel('K_bT') % x-axis label
ylabel('\langle E\rangle') % y-axis label
legend('AE20','AE40','AE60', 'AE80' )