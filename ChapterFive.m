%%  Chapter Five Examples Worked Out
% run each section seperately with [Ctrl] + [Enter] to see individual plots
clear;

%% Robustness Example
x = normrnd(1, 1, 10, 1);
% x(10) = normrnd(1, 100);
sampleMean = mean(x)
sampleVariance = var(x)
sampleMedian = median(x)

%% Example 5.27
% kernel density estimate of data shows bimodal waiting times
data = importdata('geyser.dat');
[f,xi] = ksdensity(data);
plot(xi,f,'k-','LineWidth',2)
xlabel 'Waiting Time (min)'; ylabel 'Probability Density';

%fitting mixture Gaussian distribution to this data
fun = @(lambda, data, cens, freq) - nansum(log(lambda(5)/(sqrt(2*pi)*lambda(3))*...
    exp(-((data - lambda(1))/lambda(3)).^2/2) + (1-lambda(5))/(sqrt(2*pi)*lambda(4))*...
    exp(-((data - lambda(2))/lambda(4)).^2/2)));
%nansum() is used to avoid overflow when objective function is far from its minimum
start = [55 80 5 10 .4];
parm = mle(data,'nloglf',fun,'start',start);

parm = mle(data,'nloglf',fun,'start',start,'Options',statset('MaxIter',300));
m1 = parm(1);
m2 = parm(2);
s1 = parm(3);
s2 = parm(4);
a = parm(5);

x = 20:0.1:140;
hold on
plot(x, a/(sqrt(2*pi)*s1)*exp(-((x - m1)/s1).^2/2) +...
    (1-a)/(sqrt(2*pi)*s2)*exp(-((x - m2)/s2).^2/2), 'b-','LineWidth',2);
hold off
% the Gaussian mixture model has similar modes, but is sharper than the
% kernel plot - this is probably because the GMM is not the correct model
% to use. Something with flatter peaks would be more appropriate.

%%