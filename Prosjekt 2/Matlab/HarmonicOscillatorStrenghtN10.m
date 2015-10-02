% Numerical Solution
% n = 10
data10 = load('/Users/AndreasM/ProsjekterQt/build-Project2dHarmonicOscillatorStrenghtPlot-Desktop_Qt_5_5_0_clang_64bit-Debug/ZArmadilloN10.txt');  % load eigenvector
rho10 = load('/Users/AndreasM/ProsjekterQt/build-Project2dHarmonicOscillatorStrenghtPlot-Desktop_Qt_5_5_0_clang_64bit-Debug/rhoN10.txt');  % load eigenvector

x2 = linspace(0,1,12); 


%plot(x2,v10, 'r')
%title('Numerical Solution n = 10')
%xlabel('') % x-axis label
%ylabel('') % y-axis label

title('Numerical Solution n = 10')
%plot(rho, data10(:,1))
close all
plot(rho, data10(:,1).*data10(:,1))
hold on
plot(rho, data10(:,2).*data10(:,2))
hold on
plot(rho, data10(:,3).*data10(:,3))