%% Example 9.10
% Compute the logistic regression solution on the seond annd third species,
% where y = 1 corresponds to species 2 and y = 0 corresponds to species 3

%% Importing Data
clear; close all; clf;
load fisheriris.mat

%% Data
x = [meas(51:150, :) ones(100, 1)];
y = [ones(1, 50) zeros(1, 50)]';

%% Iteratively Reweighted Logistic Regression
beta = x\y;
ovar = sum((y - x*beta).^2)/95;

for i = 1:15
    muhat = 1./(1 + exp(-x*beta));
    w = diag(muhat.*(1 - muhat));
    zeta = x*beta + (y - muhat)./(muhat.*(1 - muhat));
    beta = (x'*w*x)\(x'*w*zeta);
    nvar = sum((y - muhat).^2)/95;
    if abs(nvar - ovar) <= 0.01*ovar
        break
    end
    ovar = nvar;
end

%% Standard Errors on Coefficients
% Follows from the Fisher information matrix
sevar = x'*w*x;
se = sqrt(diag(sevar));

%% Wald Test
% Used to assess the coefficients for significance
wald = beta.^2./(se.^2);
cv = 1 - chi2cdf(wald, 1);

% The coefficents for the first three predictors accept the null hypothesis
% that the coefficients are zero. The fourth coefficient strongly rejects
% the null hypothesis.