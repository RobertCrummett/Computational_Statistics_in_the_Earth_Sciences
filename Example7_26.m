%% Example 7.26
% Analyze the combined p-values of the experiment.

%% Data
pvala = [0.04, 0.009];
pvalb = [0.22, 0.15, 0.23, 0.17, 0.20, 0.18, 0.14, 0.08, 0.31, 0.21];

%% Combined P-Values
% Two Different Methods For Data A
p_a_chi = 1 - chi2cdf(-2*sum(log(pvala)), 2*length(pvala));
p_a_norm = 1 - normcdf(sum(norminv(1 - pvala))/sqrt(length(pvala)));

% The A data rejects the null hypothesis raw and combined.

% Twho Different Methods For Data B
p_b_chi = 1 - chi2cdf(-2*sum(log(pvalb)), 2*length(pvalb));
p_b_norm = 1 - normcdf(sum(norminv(1 - pvalb))/sqrt(length(pvalb)));

% The B data does not reject the null hypothesis raw, but combined, the
% data do reject the null hypothesis.