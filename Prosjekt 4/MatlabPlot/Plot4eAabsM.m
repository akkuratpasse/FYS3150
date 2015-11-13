%Disorder
%t = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/t_Spins20.txt');
% AabsM20 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/AabsM20.txt');
% AabsM40 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/AabsM40.txt');
% AabsM60 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/AabsM60.txt');
% AabsM80 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/AabsM80.txt');
% 
% 
% plot(t, AabsM20, 'r')
% hold on
% plot(t, AabsM40, 'g')
% hold on
% plot(t, AabsM60, 'b')
% hold on
% plot(t, AabsM80, 'y')
% hold off
% title('\langle|M|\rangle: mcs = 10000, n temp = 20, Disorder')
% xlabel('K_bT') % x-axis label
% ylabel('\langle|M|\rangle') % y-axis label
% legend('AabsM20', 'AabsM40', 'AabsM60', 'AabsM80')
%Order
t = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/t.txt');
%t = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/t_Spins.txt');
AabsM20 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/AabsM20.txt');
AabsM40 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/AabsM40.txt');
AabsM60 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/AabsM60.txt');
AabsM80 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/AabsM80.txt');


plot(t, AabsM20/(20^2), 'r')
hold on
plot(t, AabsM40/(40^2), 'g')
hold on
plot(t, AabsM60/(60^2), 'b')
hold on
plot(t, AabsM80/(80^2), 'y')
hold off
title('\langle|M|\rangle: mcs = 10000, n temp = 20, ORDER')
xlabel('K_bT') % x-axis label
ylabel('\langle|M|\rangle') % y-axis label
legend('AabsM20', 'AabsM40', 'AabsM60', 'AabsM80')
