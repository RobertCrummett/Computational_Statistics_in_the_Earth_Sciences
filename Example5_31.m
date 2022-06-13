%% Example 5.31
% Find the alpha = 0.05 confidence interval on the sample median
% sizes N = 11, 101, and 1001.

%% Plotting Curves
clear; close all; clf;

Nlist = [11 101 1001]; % sample size N for each plot desired

for j = 1:length(Nlist)
    figure(j);
    hold on
    % confidence limits
    yline(0.05,'k') 
    yline(0.95,'k')

    N =Nlist(j);
    results = zeros([N 2]);
    xlim([1 N]); ylim([0 1]);

    for i = 1:N
        % Confidence intervals of the quantiles are distribution free.
        % Thompson (1936).
        results(i,:) = [i betacdf(0.5, i, N - i + 1) - betacdf(0.5, N - i + 1, 1)];  
        % Look for the conservative confidence interval.
    end

    plot(results(:,1),results(:,2),'-r')
    legend('1 - \alpha = 0.95','\alpha = 0.05',strcat("N = ",string(N)),'Location','e')
    title 'Confidence Interval of Sample Median for \alpha = 0.05'
    hold off
end