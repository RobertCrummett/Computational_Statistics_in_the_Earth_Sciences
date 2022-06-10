%% Example 4.2
% From Cavendish density of Earth data, compare the standard deviation,
% mediant absolute deviation (MAD), and interquartile range
% The truncate the data and repreat

%% Import Data
data = importdata("cavendish.dat");
n = length(data);

%% Comparison
stdBiasCorrection = sqrt(2/(n-1))*gamma(n/2)/gamma((n-1)/2);
deviations = [std(data) stdBiasCorrection*std(data,1)];
disp(deviations) %small difference between baised and unbaised std
mad(data,1) %mad larger than iqr/2 suggests symmetry
iqr(data)

%% Truncation
data = sort(data);
data = data(3:n-2); %removing 4 data points from top and bottom of sorted data

%% Comparison
1 - std(data)./deviations; 
%percent change of standard deviations when 4 data points are removed is large
mad(data,1) %mad smaller than iqr/2 suggests asymmetry increased
iqr(data)