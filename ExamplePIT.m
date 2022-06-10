%% Example: Probability Integral Transform
% Generate a set of uniform random variables and apply the inverse
% transform to get Rayleigh random variables

clear; close all; clf;

%% Sample Uniform distribution [0 1]
x = unifrnd(0,1,10000,1);

%% Apply inverse Rayleigh transform
y = raylinv(x, 1);

%% Plotting
h = histogram(y,50,'Normalization','pdf');
hold on
xx = h.BinEdges + h.BinWidth/2;
plot(xx,raylpdf(xx,1),'LineWidth',2)
hold off
xlim([0 4]); ylim([0 0.8])
title(["Probability histogram of  of 10000 draws from Raleigh distribution";...
    "using inverse prob integral transform compared with Raleigh pdf"])