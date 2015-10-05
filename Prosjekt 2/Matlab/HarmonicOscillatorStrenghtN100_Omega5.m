% Numerical Solution
% n = 100
data100_Omega5 = load('/Users/AndreasM/ProsjekterQt/build-Project2dHarmonicOscillatorStrenghtPlot-Desktop_Qt_5_5_0_clang_64bit-Debug/ZArmadilloN100_Omega5.txt');  % load eigenvector
rho100_Omega5 = load('/Users/AndreasM/ProsjekterQt/build-Project2dHarmonicOscillatorStrenghtPlot-Desktop_Qt_5_5_0_clang_64bit-Debug/rhoN100_Omega5.txt');  % load eigenvector

x2 = linspace(0,1,102); 


%plot(x2,v10, 'r')
%title('Numerical Solution n = 10')
%xlabel('') % x-axis label
%ylabel('') % y-axis label


%plot(rho, data100(:,1))
close all
plot(rho100_Omega5, data100_Omega5(:,1).*data100_Omega5(:,1))
hold on
plot(rho100_Omega5, data100_Omega5(:,2).*data100_Omega5(:,2))
hold on
plot(rho100_Omega5, data100_Omega5(:,3).*data100_Omega5(:,3))
hold on
title('Omega = 5, Numerical Solution n = 100')
xlabel('Relative coordinate, r') % x-axis label
ylabel('Radial probability, r^2') % y-axis label