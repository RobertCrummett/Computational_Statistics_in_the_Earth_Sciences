%% Example 4.11
% Obtain 10000 random draws from the Cauchy distribution, and compute a KDE
% for them. The Cauchy distribution is equivalent to the t distribution
% with one degree of freedom.

% algebraic tails on the Cauchy distribution make one expect difficulty,
% because smaller number of KDE evenly distributed acrossed the support is
% not adequate for its characterization

%% Random Draws
x = trnd(1, [1 10000]);

%% Plotting KDE with different esimators
% default number of points to evaluate is 100
[f, xi] = ksdensity(x,'npoints',100);

figure(1);
subplot(2,2,1);
plot(xi,f,'b-')
title 'npoints = 100'
subplot(2,2,2)
plot(xi,f,'b-')
title 'npoints = 100, zoomed in'
xlim([-500 500])

% now set number of points to evaluate over as 10000
[f, xi] = ksdensity(x,'npoints',10000);

subplot(2,2,3);
plot(xi,f,'r-')
title 'npoints = 10000'
subplot(2,2,4)
plot(xi,f,'r-')
title 'npoints = 10000, zoomed in'
xlim([-50 50])

% higher number of points of estimation lead to cleaner results, albeit, at
% the cost of computation time
