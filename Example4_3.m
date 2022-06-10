%% Example 4.3
% Fro Cavendish density of Earth data, compute and compare the skewness and
% kurtosis with and without applying the bias correction
% The truncate the data and repeat

%% Import data
data = importdata('cavendish.dat');

%% Compare Biased and Unbaised
%bias flag = 1 & unbias flag = 0;
skew = [skewness(data,1) skewness(data,0)];
kurt = [kurtosis(data,1) kurtosis(data,0)];
disp([skew; kurt])

%% Truncate Data
data = sort(data);
data = data(3:n-2);

%% Compare Truncated Biased and Truncated Unbiased 
skew2 = [skewness(data,1) skewness(data,0)];
kurt2 = [kurtosis(data,1) kurtosis(data,0)];
disp([skew2; kurt2])