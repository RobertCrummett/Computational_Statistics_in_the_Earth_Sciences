%% Example 11.8
% Define a sequential binary partition, and apply it to the MAR data.

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

%% Sequential Binary Partition
% Insight into the underlying scientific problem is required to suggest a
% sequential binary partition of the composition.
% Citeria used to construct the sbp are taken from Pawlosky-Glahn et al.
% (2015)
sbp = [ 0 0 0 1 -1 0 0 0 0 0 0;
    1 0 -1 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 -1;
    1 -1 1 0 0 0 0 0 0 0 -1;
    0 0 0 0 0 0 0 1 -1 0 0;
    0 0 0 0 0 0 0 1 1 -1 0;
    0 0 0 0 0 1 -1 0 0 0 0;
    0 0 0 1 1 -1 -1 0 0 0 0;
    0 0 0 1 1 1 1 -1 -1 -1 0;
    1 1 1 -1 -1 -1 -1 -1 -1 -1 1];

Psi = Contrast(sbp);
xstar = Ilr(mar, Psi);

figure(1);
boxplot(xstar)
title 'Boxplot of the Ilr computed from the MAR data'
xlabel 'Balance'

% The boxplot defines the marginal (within balance) behavior.

% Understanding the relationships between them requires the conventional
% covariance matrix and the corresponding correlation matrix.

%% Covariance and Correlation Martices
Cov = -Psi*NVariation(mar)*Psi';
ds = sqrt(diag(diag(Cov)));
Corr = inv(ds)*Cov*inv(ds);

disp(Corr)

% Examining the correlation matrix gives intuition about the between
% varaince.

%% Biplot
g = geomean(mar);
xstar = Ilr(Perturbation(mar, repmat(1./g, length(mar), 1)), Psi);
[u, s, v] = svd(xstar);
labels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};

figure(2)
subplot(1,2,1)
biplot([s(1,1)*v(:,1), s(2,2)*v(:,2)], VarLabels=labels)
axis square
title 'Covariance Biplot'

% None of the large links are close to orthogonal - therefore there is
% dependance among the balances.

subplot(1,2,2)
biplot([v(:,1), v(:,2)], VarLabels=labels)
axis square
title 'Form Biplot'

% The far from equal form biplot suggests strongly that more than two
% eigenvectors are needed to explain the data.

sgtitle 'Biplots for the Ilr Computed from the Centered MAR Data'