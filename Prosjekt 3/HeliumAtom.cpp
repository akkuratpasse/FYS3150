//This code plots the function for Project 3

x = linspace(-2,2,200); 
f=exp(-4.*sqrt(x.*x))./sqrt(x.*x);

plot(x,f, 'r')
hold on

title('Numerical Solution for Helium expectation value integral')
xlabel('0 < x < 1') % x-axis label
ylabel('v') % y-axis label
