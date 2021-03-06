% Numerical Solution
% n = 10
data10 = load('/Users/AndreasM/ProsjekterQt/build-prosjekt1bv2-Desktop_Qt_5_5_0_clang_64bit-Debug/V10.txt');
x2 = linspace(0,1,12); 

v10(2:11) = data10;
v10(1) = 0;
v10(12) = 0;

% n = 100
data100 = load('/Users/AndreasM/ProsjekterQt/build-prosjekt1bv2-Desktop_Qt_5_5_0_clang_64bit-Debug/V100.txt');
x100 = linspace(0,1,102); 

v100(2:101) = data100;
v100(1) = 0;
v100(102) = 0;

% n = 1000
data1000 = load('/Users/AndreasM/ProsjekterQt/build-prosjekt1bv2-Desktop_Qt_5_5_0_clang_64bit-Debug/V1000.txt');
x1000 = linspace(0,1,1002); 

v1000(2:1001) = data1000;
v1000(1) = 0;
v1000(1002) = 0;


% Analytical Solution
x =linspace(0,1,1000); 
exact = 1-(1-exp(-10))*x -exp(-10*x);
%plot (data10)

% Plot
%plot(x2,v,x,exact)

plot(x2,v10, 'r')
hold on
plot(x100,v100, 'g')
hold on
plot(x1000,v1000, 'b')
hold on
plot(x,exact, 'y--') % Analytical
title('Numerical and Analytical Solution for n = 10, 100, 1000')
xlabel('0 < x < 1') % x-axis label
ylabel('v') % y-axis label
