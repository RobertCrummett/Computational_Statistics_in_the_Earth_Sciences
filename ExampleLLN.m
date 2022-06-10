%% Example: Law of Large Numbers
% The mean of a large sample of random variables is close to the expected
% value of their distribution.
clear; close all; clf;

x = unidrnd(10,1000,1);
y = [];
for i = 1:1000
    y = [y mean(x(1:i))];
end

figure(1);

plot(y)
hold on
x = 5.5*ones(size(y));
plot(x,'r--')
hold off

title(["The Law of Large Numbers"; "Discrete Uniform Distribution"])
ylabel('Running Mean'); xlabel('Number of Draws');