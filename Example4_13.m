%% Example 4.13
% Compute 1000 random draws from the standard Rayleigh distribution.
% Produce varianve stabilized p-p plots using the lognormal, generalized
% extreme value, and standard Raleigh cdf.
% What are the differences, and why?

%% Random Draws
clear; close all; clf;
draws = 10005;
x = raylrnd(1, [1 draws]);
u = ((1:draws) - 0.5)/draws;
r = 2/pi*asin(sqrt(u));

%% Rayleigh Distribution
s = 2/pi*asin(sqrt(raylcdf(sort(x))));

figure(1);
subplot(1,3,1); plot(r,s,'k+')
axis square
grid on
ylabel 'Transformed Data'; xlabel 'Sine Qunatiles'; 
title 'Raleigh Distribution'

%% MLE Parameters for Lognormal Distribution
pHatLogn = lognfit(x);
s = 2/pi*asin(sqrt(logncdf(sort(x),pHatLogn(1),pHatLogn(2))));

subplot(1,3,2); plot(r,s,'k+')
axis square
grid on
xlabel 'Sine Qunatiles'; 
title 'Lognormal Distribution'

%% MLE Paramters for Generalized Extreme Value Distibution
pHatGEV = gevfit(x); %produces a [k sigma mu]
s = 2/pi*asin(sqrt(gevcdf(sort(x),pHatGEV(1),pHatGEV(2),pHatGEV(3))));

subplot(1,3,3); plot(r,s,'k+')
axis square
grid on
xlabel 'Sine Qunatiles'; 
title 'Generalized Extreme Value Distribution'

%% PDFs Plotted
%Rayleigh
figure(2);
plot(sort(x),raylpdf(sort(x),1),'-b')

%Lognormal
hold on
plot(sort(x),lognpdf(sort(x),pHatLogn(1),pHatLogn(2)),'r-')

%GEV
plot(sort(x),gevpdf(sort(x),pHatGEV(1),pHatGEV(2),pHatGEV(3)))
hold off
legend('Rayleigh PDF','Lognormal PDF','GEV PDF','Location','ne')
title 'Probability Density Functions Compared'