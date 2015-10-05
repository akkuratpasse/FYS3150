% Omega = 0.01
data100_Omega001 = load('/Users/AndreasM/ProsjekterQt/build-Project2dHarmonicOscillatorStrenghtPlot-Desktop_Qt_5_5_0_clang_64bit-Debug/ZArmadilloN100_Omega001.txt');  % load eigenvector
rho100_Omega001 = load('/Users/AndreasM/ProsjekterQt/build-Project2dHarmonicOscillatorStrenghtPlot-Desktop_Qt_5_5_0_clang_64bit-Debug/rhoN100_Omega001.txt');  % load eigenvector


% Omega = 0.5
data100_Omega05 = load('/Users/AndreasM/ProsjekterQt/build-Project2dHarmonicOscillatorStrenghtPlot-Desktop_Qt_5_5_0_clang_64bit-Debug/ZArmadilloN100_Omega05.txt');  % load eigenvector
rho100_Omega05 = load('/Users/AndreasM/ProsjekterQt/build-Project2dHarmonicOscillatorStrenghtPlot-Desktop_Qt_5_5_0_clang_64bit-Debug/rhoN100_Omega05.txt');  % load eigenvector


% Omega = 1
data100_Omega1 = load('/Users/AndreasM/ProsjekterQt/build-Project2dHarmonicOscillatorStrenghtPlot-Desktop_Qt_5_5_0_clang_64bit-Debug/ZArmadilloN100_Omega1.txt');  % load eigenvector
rho100_Omega1 = load('/Users/AndreasM/ProsjekterQt/build-Project2dHarmonicOscillatorStrenghtPlot-Desktop_Qt_5_5_0_clang_64bit-Debug/rhoN100_Omega1.txt');  % load eigenvector

% Omega = 5
data100_Omega5 = load('/Users/AndreasM/ProsjekterQt/build-Project2dHarmonicOscillatorStrenghtPlot-Desktop_Qt_5_5_0_clang_64bit-Debug/ZArmadilloN100_Omega5.txt');  % load eigenvector
rho100_Omega5 = load('/Users/AndreasM/ProsjekterQt/build-Project2dHarmonicOscillatorStrenghtPlot-Desktop_Qt_5_5_0_clang_64bit-Debug/rhoN100_Omega5.txt');  % load eigenvector


figure
% Omega = 0.01
subplot(4,1,1);
plot(rho100_Omega001, data100_Omega001(:,1).*data100_Omega001(:,1))
hold on
plot(rho100_Omega001, data100_Omega001(:,2).*data100_Omega001(:,2))
hold on
plot(rho100_Omega001, data100_Omega001(:,3).*data100_Omega001(:,3))
hold on
title('Omega = 0.01, Numerical Solution n = 100')
xlabel('r') % x-axis label
ylabel('r^2') % y-axis label


% Omega = 0.5
subplot(4,1,2);
plot(rho100_Omega05, data100_Omega05(:,1).*data100_Omega05(:,1))
hold on
plot(rho100_Omega05, data100_Omega05(:,2).*data100_Omega05(:,2))
hold on
plot(rho100_Omega05, data100_Omega05(:,3).*data100_Omega05(:,3))
hold on
title('Omega = 0.5, Numerical Solution n = 100')
xlabel('r') % x-axis label
ylabel('r^2') % y-axis label


% Omega = 1
%close all
subplot(4,1,3);
plot(rho100_Omega1, data100_Omega1(:,1).*data100_Omega1(:,1))
hold on
plot(rho100_Omega1, data100_Omega1(:,2).*data100_Omega1(:,2))
hold on
plot(rho100_Omega1, data100_Omega1(:,3).*data100_Omega1(:,3))
hold on
title('Omega = 1, Numerical Solution n = 100')
xlabel('r') % x-axis label
ylabel('r^2') % y-axis label


% Omega = 5
subplot(4,1,4);
plot(rho100_Omega5, data100_Omega5(:,1).*data100_Omega5(:,1))
hold on
plot(rho100_Omega5, data100_Omega5(:,2).*data100_Omega5(:,2))
hold on
plot(rho100_Omega5, data100_Omega5(:,3).*data100_Omega5(:,3))
hold on
title('Omega = 5, Numerical Solution n = 100')
xlabel('r') % x-axis label
ylabel('r^2') % y-axis label

