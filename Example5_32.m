%% Example 5.32
% Tolerance interval is a confidence interval for a distribution rather
% than for a statistic.

%% Estimating Tolerance Interval
% Use index i = s - r and plug into formula
clear; close all; clf;

% Inputs
N = 100;
alpha = 0.05;
plist = [0.5 0.75 0.95];
tol = zeros(length(plist),1);

for j = 1:length(plist)
    p = plist(j); % probability
    results = zeros([N 2]);

    for i = 1:N
        results(i,:) = [i (betacdf(p, N - i + 1, i) - 0.05)];
        % Close to zero implies close to boundary
    end
    [~,minIdx] = min(abs(results(:,2)));

    % Tolerance index chosen
    if results(minIdx,2) < 0
        tol(j) = results(minIdx, 1);
    else
        tol(j) = results(minIdx - 1, 1);
    end
end
