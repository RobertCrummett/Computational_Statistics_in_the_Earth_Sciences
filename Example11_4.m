%% Example 11.4
% Plot ternary diagrams of the oxides (Fe2O3, K2O, P2O5) that exhibit the
% highest covariation in Example 11.3 before and after centering

%% Importing Data
clear; clf; close all;
mar = readmatrix("MAR.xls");

%% Data Cleaning
% Keep only the major oxides
mar = mar(:, 3:13);

% Keep only data will all of the major oxides
mar = [mar(1:4, :)', mar(22, :)', mar(47:88, :)']';

% Remove redundant outliers at 37-38 and 39-40
% Remove redundant data at 24-26 and 30-36
mar = [mar(1:25, :)', mar(27:35, :)', mar(41:47, :)']';

mar = Close(mar, 100);

%% Ternary Diagrams Before and After Centering
% Before Centering
mar1 = Close(mar(:,[4, 10, 11]), 100);

figure(1)
subplot(1,2,1)
Ternary;
TernaryPlot(mar1, 'd', 'k')
TernaryLabel('Fe2O3','K2O','P2O5')

% After Centering
subplot(1,2,2)
g = Close(geomean(mar1));
Ternary(0, 1./g, [0 0.1 0.5 0.9 1])
cmar1 = Perturbation(mar1, repmat(1./g, length(mar1), 1));
TernaryPlot(cmar1, 'd', 'b')
TernaryLabel('Fe2O3','K2O','P2O5')

sgtitle('Ternary Plots Before and After Centering Data')
