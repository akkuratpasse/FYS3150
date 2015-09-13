% Numerical Solution
data = load('/Users/AndreasM/ProsjekterQt/build-prosjekt1bv2-Desktop_Qt_5_5_0_clang_64bit-Debug/V.txt');
x2 = linspace(0,1,12); 

v(2:11) = data;
v(1) = 0;
v(12) = 0;

% Analytical Solution
x =linspace(0,1,1000); 
exact = 1-(1-exp(-10))*x -exp(-10*x);
plot (data)

% Plot
%plot(x2,v,x,exact)

plot(x2,v, 'r')
hold on
plot(x,exact, 'b')
