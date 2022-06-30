%% Example 11.9
% Komatiites are ultramafic mantle-derived rocks of predominantly Archaean
% age that exhibit compositions that are consistent with melting at
% temperatures of over 1600 degrees Celcius. Lets explore these raskals.

%% Importing Data
clear; clf; close all;
x = importdata("komatiite.dat");

%% Data
komatiite  = [x(:, 1:4) x(:, 10) x(:, 5:6) x(:, 9) x(:, 8) x(:, 7) x(:, 12) x(:, 11) x(:, 13)];
komatiite = Close(komatiite);

center = Close(geomean(komatiite), 100);

[T, totvar] = NVariation(komatiite);

%% Eigenvalues Calculated
g = Close(geomean(komatiite));
[u, s, v] = svd(Clr(Perturbation(komatiite, repmat(1./g, length(komatiite), 1))));

% Eigenvalues
disp(diag(s).^2./sum(diag(s).^2)')

%Scree Plot
figure(1)
plot(diag(s).^2./sum(diag(s).^2), 'r-', 'LineWidth', 1.5)
xlabel 'Index'; ylabel 'Normalized Eigenvalue'
title 'Scree Plot for the SVD of the Clr of the Centered Komatiite Data'

% This suggests that at least five eigenvalues are needed to explain the
% data.

disp(v(:,1:9))

% The mantle signature of Cr vs Ni only appears in the fifth or sixth
% eigenvector. This suggests that metamorphism in the past two billion
% years is responsible for most of the variation, and that these clr
% eigenvectors will not provide gerat insight.

%% Sequential Binary Partition

sbp = [1 0 -1 0 0 0 0 0 0 0 0 0 0;
    0 1 0 0 -1 0 0 0 0 0 0 0 0;
    0 1 0 0 1 0 0 0 0 0 0 0 -1;
    1 -1 1 0 -1 0 0 0 0 0 0 0 -1;
    0 0 0 0 0 0 0 0 0 1 -1 0 0;
    0 0 0 0 0 0 0 0 0 1 1 -1 0;
    0 0 0 0 0 1 -1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 -1 0 0 0 0;
    0 0 0 0 0 1 1 -1 -1 0 0 0 0;
    0 0 0 1 0 -1 -1 -1 -1 0 0 0 0;
    0 0 0 1 0 1 1 1 1 -1 -1 -1 0;
    1 1 1 -1 1 -1 -1 -1 -1 -1 -1 -1 1];

Psi = Contrast(sbp);
xstar = Ilr(komatiite, Psi);

figure(2)
boxplot(xstar)
xlabel 'Balance'
title 'Boxplot of the Ilr of the Komatiite Data'

%% Covariance and Correlation Matrix
Cov = -Psi*NVariation(komatiite)*Psi';
ds = sqrt(diag(diag(Cov)));
Corr = inv(ds)*Cov*inv(ds);

%% Biplots
g = geomean(komatiite);
xstar = Ilr(Perturbation(komatiite, repmat(1./g, length(komatiite), 1)), Psi);
[u, s, v] = svd(xstar);

labels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'};
figure(3)
biplot([s(1,1)*v(:,1), s(2,2)*v(:,2)], 'VarLabels', labels)
title 'Covariance Biplot of the Ilr of the Centered Komatiite Data'

figure(4)
biplot([s(2,2)*v(:,2), s(3,3)*v(:,3)], 'VarLabels', labels)
title 'Covariance Biplot of the Ilr of the Centered Komatiite Data'
xtitle 'Component 2'; ylabel 'Component 3';