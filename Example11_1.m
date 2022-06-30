%% Example 11.1
% Using Mid-Ocean Ridge (MAR) data, demonstrate their inconsistenci with
% standard estimators for correction and the requirement to use
% compositional concepts.

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

%% Examine Correlations
% In classical statistics, the standard descriptive statistics are the
% sample mean and variance (or covaraince). However, these are not suitable
% measures of central tendancy and dispersion for compositions, as can be
% readily shown.
% Respecively, the correlation coefficents displayed are (CaO, K2O),
% (Fe2O3, FeO), (FeO, MgO), (Fe2O3, MgO)
disp([corr(mar(:, 8), mar(:, 10)); corr(mar(:, 4), mar(:, 5)); corr(mar(:, 5), mar(:, 7)); corr(mar(:, 4), mar(:, 7))]);

% Delete SiO2 and Al2O3 from the data set, reclose and recompute the
% correlation coefficients.

mar1 = [mar(:,2) mar(:, 4:11)];
mar1 = Close(mar1, 100);
disp([corr(mar1(:, 6), mar1(:, 8)); corr(mar1(:, 2), mar1(:, 3)); corr(mar1(:, 3), mar1(:, 5)); corr(mar1(:, 2), mar1(:, 5))]);

% The correlation matrix between and mafic constituents (MgO, FeO, Fe2O3)
% before closing the subcomposition is:

disp(corr([mar(:, 4:5) mar(:, 7)]))

% And after closing the sub composition is:
mar1 = Close([mar(:, 4:5) mar(:, 7)], 100);
disp(corr(mar1))

% The differences are dramatic and illustrate the incompatibility of
% classical summary statistics with the Aitchison geometry.

