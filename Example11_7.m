%% Example 11.7
% Compute covariance and form biplots for the first two eigenvectors of the
% MAR data, and present the key insights they provide as ternary diagrams
% before and after centering.

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

%% Singular Value Decomposition
g = Close(geomean(mar));
[u, s, v] = svd(Clr(Perturbation(mar, repmat(1./g, length(mar), 1))));

%% Biplots
labels = {'SiO2','TiO2','Al2O3','Fe2O3','FeO','MnO','MgO','CaO','Na2O','K2O','P2O5'};

figure(1);
subplot(1,2,1)
biplot([s(1,1)*v(:,1), s(2,2)*v(:,2)], 'VarLabels', labels);
axis square
title 'Covariance Biplot'

% This biplot is somewhat one dimensional, reflecting the size of the first
% singular value.
% The rays for K2O, Fe2O3, and P2O5 dominate.
% The rays for SiO2, Al2O3, MnO, CaO are nearly coincident, and ar thus
% redundant. They coincide in direction but not magnitude for FeO and MgO.
% The links between FeO-P2O5 and Fe2O3-P2O5 are nearly orthogonal, so their
% ratios are independant. 

subplot(1,2,2)
biplot([v(:,1), v(:,2)], 'VarLabels', labels);
axis square
title 'Form Biplot'

sgtitle 'Biplot for the MAR Oxide Data'

% The length of the predominant oxides are approximately equal, but the
% entire data set is not evenly represented. This suggests that more
% eigenvalues may be needed.

figure(2);
subplot(1,2,1)
biplot([s(1,1)*v(:,1), s(2,2)*v(:,2), s(3,3)*v(:,3)], 'VarLabels', labels);
axis square
title 'Covariance Biplot'

subplot(1,2,2)
biplot([v(:,1), v(:,2), v(:,3)], 'VarLabels', labels);
axis square
title 'Form Biplot'

sgtitle 'Biplot for the MAR Oxide Data'

% The form biplot is more evenly distributed with three singular vectors,
% however, the covaraince biplot is really hard to evaluate because so many
% components are nearly coincident.

%% Ternary Diagrams
mar1 = Close(mar(:,[5, 9, 11]), 100);

figure(3)
subplot(1,2,1)
Ternary;
TernaryPlot(mar1, 'd', 'k');
TernaryLabel(labels{5},labels{9},labels{11})

subplot(1,2,2)
g = Close(geomean(mar1));
Ternary(0, 1./g, [0 0.1 0.3 0.7 0.9]);
cmar1 = Perturbation(mar1, repmat(1./g, length(mar1), 1));
TernaryPlot(cmar1, 'd', 'b');
TernaryLabel(labels{5},labels{9},labels{11})

sgtitle 'Ternary Diagrams for FeO - Na_{2}O - P_{2}O_{5} Before and After Centering'

% The ternary diagram shows little correlation before or after centering,
% consistent with what what seen in the biplots - these components had
% links at approximately right angles.

mar2 = Close(mar(:,[4, 5, 7]), 100);

figure(4)
subplot(1,2,1)
Ternary;
TernaryPlot(mar2, 'd', 'k');
TernaryLabel(labels{4},labels{5},labels{7})

subplot(1,2,2)
g = Close(geomean(mar2));
Ternary(0, 1./g, [0 0.1 0.3 0.5 0.7 0.9]);
cmar1 = Perturbation(mar2, repmat(1./g, length(mar2), 1));
TernaryPlot(cmar1, 'd', 'b');
TernaryLabel(labels{4},labels{5},labels{7})

sgtitle 'Ternary Diagram for the Ferric Minerals Before and After Centering'

% The biplots suggess that there is a compositional line between these
% components. This is really obvious in the centered data.

mar3 = Close(mar(:,[8, 9, 10]), 100);

figure(5)
subplot(1,2,1)
Ternary;
TernaryPlot(mar3, 'd', 'k');
TernaryLabel(labels{8},labels{9},labels{10})

subplot(1,2,2)
g = Close(geomean(mar3));
Ternary(0, 1./g, [0 0.2 0.8 1]);
cmar1 = Perturbation(mar3, repmat(1./g, length(mar3), 1));
TernaryPlot(cmar1, 'd', 'b');
TernaryLabel(labels{8},labels{9},labels{10})

sgtitle 'Ternary Diagram of the Feldspar Minerals Before and After Centering'

% The centered data shows that there are nearly constant proportions of 0.8
% CaO to 0.2 Na2O and K2O