data = load('/Users/AndreasM/ProsjekterQt/build-prosjekt1bv2-Desktop_Qt_5_5_0_clang_64bit-Debug/V.txt');
x2 = linspace(0,1,12);
x =linspace(0,1,1000); 
exact = 1-(1-exp(-10))*x -exp(-10*x)  % exact
plot (data)

%x = 0:9;
v(2:11) = data;
v(1) = 0;
v(12) = 0;
plot(x2,v,x,exact)
y = V;
u = U;
% hold on
% plot (x,y, 'Color',[0,0.7,0.9]);
% plot (x,u);
% title('Solution n = 10')
% xlabel('x')
% ylabel('v(x)')

