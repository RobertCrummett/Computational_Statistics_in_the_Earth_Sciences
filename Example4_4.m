%% Example 4.4
% The file paleocurrent.dat contains 30 values of the cross-bed azimuth of
% paleocurrents measued in the Bedford Anticline of Australia. Compute the
% sample mean direction, mean resultant length, circular variance, smaple
% skewness, and the sample kurtosis.
% What do these tell us about the data? Fit a Von Mises distribution to the
% data and qualitatively evaluate the fit.

%% Import Data
clear; clf;
data = importdata("paleocurrent.dat");
figure(1);
polarhistogram(pi*data/180, 'NumBins',10) % clusted around $250^{\circ}$
title 'Rose diagram of the paleocurrent data'

%% Compute Vector Resultants
chat = sum(cosd(data));
shat = sum(sind(data));
thetabar = atand(shat/chat); %sample mean direction
thetabar = thetabar + 180; %shat and chat indicate thetabar is in third quadrant

%% Resultant Length
rhat = sqrt(chat^2 + shat^2);
rbar = rhat/length(data);
vhat = 1 - rbar; %sample circular variance

%% kth Central Trignometric Moment
% k =2
c2hat = sum(cosd(2*(data - thetabar)));
s2hat = sum(sind(2*(data - thetabar)));

theta2bar = atand(s2hat/c2hat);
r2bar = sqrt(c2hat^2 + s2hat^2)/length(data);

skew = r2bar*sind((theta2bar - 2*thetabar))/(1 - rbar)^1.5; 
% sample skewness for circular data
kurt = (r2bar*cosd(theta2bar - 2*thetabar) - rbar^4)/(1 - rbar)^2; 
% sample kurtosis for circualr data

% skewness and kurtosis suggest a dist. that is heavier toward samller
% angles and platykurtic

%% Von Mises distribution fit
% takes direction theta bar and the concentration parameter can be found
% with the non linear root finder fzero()

fun = @(x) besseli(1, x)./besseli(0, x) - rbar;
% solve for where this equals zero

kappabar = fzero(fun,1);
% the maximum likelihood estimator of the concentration parameter

%% Ploting Von Mises distribution
figure(2);
histogram(data,10,'Normalization','pdf')

hold on
theta = 160:0.1:340;
plot(theta,exp(kappabar*cosd(theta - thetabar))/(360*besseli(0,kappabar)),'LineWidth',2)
xlim([160 330])
ylabel('Probablility Density'); xlabel('Angle (degrees)');
title 'Empiracle pdf of paleocurrent data compared to a Von Mises distribution'
hold off
%empirical pdf of paleocurrent data compared to a Von Mises distribution
%with parameters kappabar = 0.26757, mubar = 247.6166

