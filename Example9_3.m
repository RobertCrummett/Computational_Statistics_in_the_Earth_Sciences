%% Example 9.3
% Repeat the analysis of the star data using the bootstrap to estimate the
% coefficents and test them for significance.

%% Importing Data
clear; close all; clf;
star = importdata('star.dat');

%% Data
y = star(:, 2);
x = [star(:, 1) ones(size(y))];

%% Bootstrapping
rng default
boot = bootstrp(99999, @(x, y) x\y, x, y);

params = mean(boot);

% The bootstrap distributions are skewed by the leverage points, with a
% long upper tail for the slope and a long lower tail for the intercept.

%% Plot Bootstrap Distributions
figure(1);
subplot(1,2,1)
[f, xi] = ksdensity(boot(:,1));
plot(xi, f, 'b-','LineWidth', 2)
title 'Slope Coefficient'

subplot(1,2,2)
[f, xi] = ksdensity(boot(:,2));
plot(xi, f, 'b-', 'LineWidth', 2)
title 'Intercept Coefficient'

sgtitle 'Bootstrap Distributions for the Coefficients of the Star Data'

%% Implement the Double Bootstrap
b = x\y;
r = y - x*b;
sd = sqrt(r'*r/(length(y) - 2)*diag(inv(x'*x)));
that = b'./sd'; % test statistic
boot = bootstrp(999, @(x, y) func(x, y, b'), x, y);

for i = 1:2
    % for the significance of the fit
    pval2(i) = 2*min(sum(boot(:,i) > that(i)), sum(boot(:,i) <= that(i)))/999; 
end

% The slope parameter is not significant, but the intercept is required.
% Delete the four leverage points and repeat.

%% Plot Bootstrap Distributions for Assessing Significance of Fit
figure(2);
subplot(1,2,1)
[f, xi] = ksdensity(boot(:,1));
plot(xi, f, 'r-','LineWidth', 2)
hold on
xline(that(1),'k--')
hold off
title 'Slope'
xlim([-5 10])

subplot(1,2,2)
[f, xi] = ksdensity(boot(:,2));
plot(xi, f, 'r-', 'LineWidth', 2)
hold on
xline(that(2),'k--')
hold off
title 'Intercept'
xlim([-10 10])

sgtitle 'Bootstrap Distributions for the Significance of Fit'

%% Truncating the Star Data
y([7 11 20 30 34 36]') = [];
x([7 11 20 30 34 36]', :) = [];

b = x\y;
r = y - x*b;

sd = sqrt(r'*r/(length(y) - 2) * diag(inv(x'*x)));
that2 = b'./sd';

boot = bootstrp(9999, @(x,y) (x'*x)\x'*y, x, y);
params2 = mean(boot);

for i = 1:2
    pval3(i) = 2*min(sum(boot(:,i) > that2(i)), sum(boot(:,i) <= that2(i)))/999; 
end


%% Plotting Truncated Bootstrap Distributions
figure(3);
subplot(1,2,1)
[f, xi] = ksdensity(boot(:,1));
plot(xi, f, 'b-','LineWidth', 2)
title 'Slope Coefficient'

subplot(1,2,2)
[f, xi] = ksdensity(boot(:,2));
plot(xi, f, 'b-', 'LineWidth', 2)
title 'Intercept Coefficient'

sgtitle 'Bootstrap Distributions for the Coefficients of the Truncated Star Data'

%% Implement the Double Bootstrap to Truncated Data
sd = sqrt(r'*r/(length(y) - 2)*diag(inv(x'*x)));
boot = bootstrp(999, @(x, y) func(x, y, b'), x, y);

for i = 1:2
    % for the significance of the fit
    pval4(i) = 2*min(sum(boot(:,i) > that2(i)), sum(boot(:,i) <= that2(i)))/999; 
end

%% Plot Bootstrap Distributions for Truncated Data for Assessing Significance of Fit
figure(4);
subplot(1,2,1)
[f, xi] = ksdensity(boot(:,1));
plot(xi, f, 'r-','LineWidth', 2)
hold on
xline(that2(1),'k--')
hold off
title 'Slope'
xlim([-10 10])

subplot(1,2,2)
[f, xi] = ksdensity(boot(:,2));
plot(xi, f, 'r-', 'LineWidth', 2)
hold on
xline(that2(2),'k--')
hold off
title 'Intercept'
xlim([-10 10])

sgtitle 'Bootstrap Distributions for Truncated Data for the Significance of Fit'

%% Function
function [t] = func(x, y, b)
opt = statset('UseParallel',true); % Parallelize the bootstrap.
bmat = bootstrp(999, @(x, y) x\y, x, y, 'Options', opt);
sd = std(bmat, 1);
t = (mean(bmat) - b)./sd;
end