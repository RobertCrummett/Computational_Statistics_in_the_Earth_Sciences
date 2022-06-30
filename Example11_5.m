%% Example 11.5
% Present and describe a dendrogram for the clr of the centered MAR data.

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

%% Dendrogram
labels = {'SiO2','TiO2','Al2O3','Fe2O3','FeO','MnO','MgO','CaO','Na2O','K2O','P2O5'};
g = Close(geomean(mar));
tree = linkage(Clr(Perturbation(mar, repmat(1./g, length(mar), 1)))');
dt = dendrogram(tree, 'Labels', labels);
title 'Dendrogram of the Clr of the Centered MAR Data'

% It is interesting that the incompatible elements Fe3+, K, and P cluster
% at the right and are together a large distance from the remaining ones
% but with Ti the closest of them,  followed by Mn. The remaining elements
% are close together in the clusters (Fe2+, Mg, Na) and (Si, Al, Ca).