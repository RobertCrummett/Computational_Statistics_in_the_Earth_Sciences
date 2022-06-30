%% Figure 11.1
clear; clf; close all;
% First Set of Parallel Lines
x = [(0:0.01:1)' (0:0.01:1)' 0.5*ones(101, 1)];
Ternary;

x0 = repmat([1/3 1/3 1/3], size(x, 1), 1);
alpha = 10;
Result = Perturbation(x0, Powering(alpha, x));
TernaryLine(Result, 'b')

x0 = repmat([1/6 1/2 1/3], size(x, 1), 1);
Result = Perturbation(x0, Powering(alpha, x));
TernaryLine(Result, 'b')

x0 = repmat([1/2 1/6 1/3], size(x, 1), 1);
Result = Perturbation(x0, Powering(alpha, x));
TernaryLine(Result, 'b')

% Second Set of Parallel Lines
x = [(0:0.01:1)' 0.5*ones(101, 1) (1:-0.01:0)'];

x0 = repmat([1/3 1/3 1/3], size(x, 1), 1);
Result = Perturbation(x0, Powering(alpha, x));
TernaryLine(Result, 'r')

x0 = repmat([1/6 1/2 1/3], size(x, 1), 1);
Result = Perturbation(x0, Powering(alpha, x));
TernaryLine(Result, 'r')

x0 = repmat([1/2 1/6 1/3], size(x, 1), 1);
Result = Perturbation(x0, Powering(alpha, x));
TernaryLine(Result, 'r')