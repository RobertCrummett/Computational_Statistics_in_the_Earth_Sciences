%% Example 4.15
% Compute 1000 random draws from the standard Rayleigh distribution.
% Produce q-q plots using the lognormal and generalized extreme value (GEV)
% distributions with maximum likelihood estimator (MLE) parameters and the
% standard Rayleigh cdf. What differences do you see, and why?

%% Random Draws
clear; close all; clf;
n = 1000; % number of draws
x = raylrnd(1,1,n);
u = ((1:n) - 0.5)/n;

%% Plotting q-q plots
% Standard Rayleigh quanitles
figure(1);
subplot(1,3,3); plot(raylinv(u,1),sort(x),'k+')
axis square
grid on
ylim([0 4]);
xlabel 'Rayleigh Qunatiles';

% MLE Paramters for GEV
pHat = gevfit(u);
% GEV quantiles
subplot(1,3,2); plot(gevinv(u,pHat(1),pHat(2),pHat(3)),sort(x),'k+')
axis square
grid on
ylim([0 4]); xlim([-5 5]);
xlabel 'GEV Qunatiles';

% MLE Parameters for Lognormal
pHat = lognfit(u);
% Lognormal quanitles
subplot(1,3,1); plot(logninv(u,pHat(1),pHat(2)),sort(x),'k+')
axis square
grid on
ylim([0 4]);
xlabel 'Log Normal Qunatiles'; ylabel 'Sorted Data';

% The lognormal distribution initially has a surplus of probability,
% followed by a deficit. This is clear from the inital steep positive slope
% and then decrease.
% The GEV distiribution exists for negative arguments, which is masked in
% the p-p plots. An initial deficit of probability exists, as is evidencesd
% by the initial upward curvature.
% Rayleigh is a straight line.