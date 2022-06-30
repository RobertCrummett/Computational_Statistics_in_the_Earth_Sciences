%% Example 10.2
% Zonally averaged temperature data for various ranges and latitiudes from
% the 1880's to the present are available at
% http://data.giss.nasa.gov/gistemp/
% The data are available with and without sea surface temperature (SST) and
% are presented as a temperature anomaly relative to a base period of
% 1951 - 1980. The data that include SST will be evaluated for significance
% of change over time by binning them into age groups and rpeforming a
% MANOVA test.

%% Importing Data
clear; clf; close all;
x = importdata('ZonAnnTs+dSST.dat'); 

% N and S hemisphere data in columns 3 and 4

%% Data
% Plotting Data
figure(1);
plot(1880:2015, x(:, 3), 'b-', 'LineWidth', 1.25)
hold on
plot(1880:2015, x(:, 4), 'r-', 'LineWidth', 1.25)
hold off

ylim([-100 150])
xlabel 'Year'
ylabel 'Temperature  Anomaly'
title 'Zonally Averaged Temperature Anomaly'
legend 'Northern Hemisphere' 'Southern Hemisphere' 'Location' 'se'

p = 2; % categroies (N and S) 
m = 13; % number of bins
n = 10; % number of samples per bin
xh = zeros(n, m, p); % preallocation
i = 1:n;

% Rearranging into equal bin sizes
for j = 1:m
    xh(i, j, 1:p) = x((j - 1)*n + i, 3:4);
end

% Preallocation
e = zeros(p,p);
h = zeros(p,p);

%% Wilk's Lambda Distribution
% Hypothesis (between) and Error (within) Matrices
for i = 1:p
    for j = 1:p
        h(i, j) = sum(sum(xh(:, :, i)*(eye(m) - ones(m, m)/m)*xh(:, :, j)'))/n;
        e(i, j) = trace(xh(:, :, i)*xh(:, :, j)') - sum(sum(xh(:, :, i)*xh(:, :, j)'))/n;
    end
end

% Wilk's Lambda Distribution
[pval, lambda] = WilksLambda(h, e, p, m - 1, n*(m - 1));

% Since the p-value is zero, the null hypothesis can be rejected that the
% means of N and S hemisphere data are the same over time.

%% Plot Latitude Zonation
% Northern Zonation
figure(2)
ndata = 8:11;
clr = {'m' 'b' 'g' 'k'};

hold on
for i = 1:length(ndata)
    plot(1880:2015, x(:,ndata(i)),'Color', clr{i}, 'LineWidth', 1.25)
end
hold off

ylim([-200 300])
xlabel 'Year'
ylabel 'Temperature  Anomaly'
title 'Zonally Averaged Temperature Anomaly in Northern Hemisphere'
legend 'Latitude: 0 - 24' 'Latitude: 24 - 44' 'Latitude: 44 - 64' 'Latitude: 64 - 90' Location northeastoutside

% Southern Zonation
figure(3)
ndata = 12:15;
clr = {'m' 'b' 'g' 'k'};

hold on
for i = 1:length(ndata)
    plot(1880:2015, x(:,ndata(i)),'Color', clr{i}, 'LineWidth', 1.25)
end
hold off

ylim([-300 200])
xlabel 'Year'
ylabel 'Temperature  Anomaly'
title 'Zonally Averaged Temperature Anomaly in Southern Hemisphere'
legend 'Latitude: 0 - 24' 'Latitude: 24 - 44' 'Latitude: 44 - 64' 'Latitude: 64 - 90' Location northeastoutside

% There is much higher variablility at the polar latitudes.

%% Wilk's Lambda Distribution for Northern/Southern Zontations
decision = 3; % 2 for both North and South, 1 for North, 0 for South, 3 to exclude poles

m = 13; % number of bins
n = 10; % number of samples per bin
i = 1:n;

if decision == 2
    p = 8;
    xh = zeros(n, m, p); % preallocation
    % Northern and Southern Data
    for j = 1:m
        xh(i, j, 1:p) = x((j - 1)*n + i, 8:15);
    end
elseif decision == 1
    p = 4;
    xh = zeros(n, m, p); % preallocation
    % Northern Data
    for j = 1:m
        xh(i, j, 1:p) = x((j - 1)*n + i, 8:11);
    end
elseif decision == 0
    p = 4;
    xh = zeros(n, m, p); % preallocation
    % Southern Data
    for j = 1:m
        xh(i, j, 1:p) = x((j - 1)*n + i, 12:15);
    end
elseif decision == 3
    p = 6;
    xh = zeros(n, m, p); % preallocation
    % Exlucing Poles from the Data
    for j = 1:m
        xh(i, j, 1:p) = x((j - 1)*n + i, 9:14);
    end
end

% Preallocation
e = zeros(p,p);
h = zeros(p,p);

% Hypothesis (between) and Error (within) Matrices
for i = 1:p
    for j = 1:p
        h(i, j) = sum(sum(xh(:, :, i)*(eye(m) - ones(m, m)/m)*xh(:, :, j)'))/n;
        e(i, j) = trace(xh(:, :, i)*xh(:, :, j)') - sum(sum(xh(:, :, i)*xh(:, :, j)'))/n;
    end
end

% Wilk's Lambda Distribution
[pval, lambda] = WilksLambda(h, e, p, m - 1, n*(m - 1));

% For the N hemisphere data, the null hypothesis that the 10 year averages
% are the same is rejected.

% For the N & S hemisphere data, the p-value is again 0, resulting in a
% rejection of the null hypothesis that the 10 year averages are the same.

% Thus the data set shows that mean temperature has changed over time.
% However, the change is far from just an increase. There are spacial
% patterns that cannot be discerned by just observing the means.
