
% ORDER
OT1 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/ORDERenergiesT=1.000000.txt');
OT24 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/ORDERenergiesT=2.400000.txt');
% DISORDER
DT1 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/DISORDERenergiesT=1.000000.txt');
DT24 = load('/Users/AndreasM/ProsjekterQt/build-Prosjekt4d-Desktop_Qt_5_5_0_clang_64bit-Release/DISORDERenergiesT=2.400000.txt');

figure
subplot(2,2,1)
histogram(OT1,5);
xlabel('\langle E\rangle')
title('OT1 = Order')
%legend('OT1 = Order')

subplot(2,2,2)
histogram(DT1,5);
xlabel('\langle E\rangle')
title('DT1 = Disorder')
%legend('DT1 = Disorder')

subplot(2,2,3)
histogram(OT24);
xlabel('\langle E\rangle')
title('OT24 = Order')
%legend('OT24 = Order')

subplot(2,2,4)
histogram(DT24);
xlabel('\langle E\rangle')
title('DT24 = Disorder')
%legend('DT24 = Disorder')


[a, b] = hist(rand(5000,1),50);
bar(a,b)
XData cannot contain duplicate values.
 
bar(b,a)
close all
bar(b,a)
bar(b/sum(b),a)
close all
bar(b/sum(b),a)
bar(b,a/sum(a))
close all
bar(b,a/sum(a))
Plot4cTempOrderDisorder