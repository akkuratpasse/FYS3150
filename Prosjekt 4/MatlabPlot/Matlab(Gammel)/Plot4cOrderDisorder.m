% Mcs = 10 000

% High T
% mcs = 10000
x = linspace(0,1,10000);
AE = load('/Users/AndreasM/ProsjekterQt/build-Project4cOrdervsDisorder-Desktop_Qt_5_5_0_clang_64bit-Debug/AE.txt');
E1 = load('/Users/AndreasM/ProsjekterQt/build-Project4cOrdervsDisorder-Desktop_Qt_5_5_0_clang_64bit-Release/energiesT=1.000000.txt');

plot(x, AE, 'r')
title('High T = 2.4 mcs= 10')
xlabel('') % x-axis label
ylabel('') % y-axis label

% Low T

