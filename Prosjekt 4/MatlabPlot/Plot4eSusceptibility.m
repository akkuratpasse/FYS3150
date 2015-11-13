% t = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/t_Spins20.txt');
% 
% %Disorder
% Susceptibility20 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/Susceptibility20.txt');
% Susceptibility40 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/Susceptibility40.txt');
% Susceptibility60 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/Susceptibility60.txt');
% Susceptibility80 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4e-Desktop_Qt_5_5_0_clang_64bit-Release/Susceptibility80.txt');
% 
% plot(t, Susceptibility20, 'r')
% hold on
% plot(t, Susceptibility40, 'g')
% hold on
% plot(t, Susceptibility60, 'b')
% hold on
% plot(t, Susceptibility80, 'y')
% hold off
% title('Susceptibility: mcs = 10000, n temp = 20, Disorder')
% xlabel('k_bT') % x-axis label
% ylabel('Susceptibility') % y-axis label
% legend('Susceptibility20','Susceptibility40', 'Susceptibility60','Susceptibility80'  )

%Order
t = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/t.txt');
Susceptibility20 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/Susceptibility20.txt');
Susceptibility40 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/Susceptibility40.txt');
Susceptibility60 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/Susceptibility60.txt');
Susceptibility80 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4bV5-Desktop_Qt_5_5_0_clang_64bit-Release/Susceptibility80.txt');


plot(t, Susceptibility20, 'r')
hold on
plot(t, Susceptibility40, 'g')
hold on
plot(t, Susceptibility60, 'b')
hold on
plot(t, Susceptibility80, 'y')
hold off
title('Susceptibility: mcs = 10000, ntemp = 20, ORDER')
xlabel('k_bT') % x-axis label
ylabel('Susceptibility') % y-axis label
legend('Susceptibility20','Susceptibility40', 'Susceptibility60','Susceptibility80'  )

