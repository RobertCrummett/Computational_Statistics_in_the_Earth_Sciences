%% Example: F Distributions
% Arises where the variance of two independant variables are compared.
% There is a tendancy toward Gaussianity the degrees of freedom for
% both random variables get large. 
% The order of parameters is important. 

%% Central F PDF's compared
% mu and nu are the degrees of freedom for the chi squared distributed
% random variables
clear; close all; clf;
x = 0:0.01:4;

figure(1);
plot(x, fpdf(x,1,1),'b')
xlim([0 4]); ylim([0 4]);
ylabel 'F PDF'
title 'F PDF for paired degrees of freedom'

hold on
plot(x, fpdf(x,2,10),'r')
plot(x,fpdf(x,10,100),'g')
plot(x, fpdf(x,100,100),'m')
hold off

legend('\mu = 1, \nu = 1','\mu = 2, \nu = 10','\mu = 10, \nu = 100',...
    '\mu = 100, \nu = 100','Location','ne')

%% Non central F distributions compared to Central F PDF
% lambda is the non centrality parameter (MATLAB refers to this paramter as
% delta in their functions)
x = 0:0.01:5;

figure(2);
plot(x, fpdf(x,10,100),'b')

hold on
plot(x, ncfpdf(x,10,100,1),'r')
plot(x,ncfpdf(x,10,100,3),'g')
plot(x,ncfpdf(x,10,100,10),'m')
plot(x,ncfpdf(x,10,100,30),'c')
hold off

ylabel 'Probability Density'
title 'Central F distribution (\mu = 10, \nu = 100) compared to Noncentral F distributions (\mu = 10, \nu = 100)'
legend('Central F PDF','\lambda = 1','\lambda = 3','\lambda = 10',...
    '\lambda = 30','Location','ne')
