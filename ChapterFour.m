%%  Chapter Four Examples Worked Out
% run each section seperately with [Ctrl] + [Enter] to see individual plots
clear;

%% Example 4.1
data = importdata('cavendish.dat');
histogram(data,10)
means = [mean(data);
harmmean(data);
geomean(data);
trimmean(data, 5)];
disp(means)

%% Example 4.2
n = length(data);
stdBiasCorrection = sqrt(2/(n-1))*gamma(n/2)/gamma((n-1)/2);
deviations = [std(data) stdBiasCorrection*std(data,1)];
disp(deviations) %small difference between baised and unbaised std
mad(data,1) %mad larger than iqr/2 suggests symmetry
iqr(data)

data = sort(data);
data = data(3:n-2); %removing 4 data points
1 - std(data)./deviations;
%percent change of standard deviations when 4 data points are removed is large
mad(data,1) %mad smaller than iqr/2 suggests asymmetry increased
iqr(data)

%% Example 4.3
data = importdata('cavendish.dat');
%for everything, bias flag = 1 & unbias flag = 0;
skew = [skewness(data,1) skewness(data,0)];
kurt = [kurtosis(data,1) kurtosis(data,0)];
disp([skew; kurt])

data = sort(data);
data = data(3:n-2);
skew2 = [skewness(data,1) skewness(data,0)];
kurt2 = [kurtosis(data,1) kurtosis(data,0)];
disp([skew2; kurt2])

%% Example 4.4
data = importdata("paleocurrent.dat");
rose(pi*data/180)

chat = sum(cosd(data));
shat = sum(sind(data));
thetabar = atand(shat/chat);
thetabar = thetabar + 180; %shat and chat indicate thetabar is in third quadrant

rhat = sqrt(chat^2 + shat^2);
rbar = rhat/length(data);
vhat = 1 - rbar;
c2hat = sum(cosd(2*(data - thetabar)));
s2hat = sum(sind(2*(data - thetabar)));
theta2bar = atand(s2hat/c2hat);
r2bar = sqrt(c2hat^2 + s2hat^2)/length(data);
skew = r2bar*sind((theta2bar - 2*thetabar))/(1 - rbar)^1.5;
kurt = (r2bar*cosd(theta2bar - 2*thetabar) - rbar^4)/(1 - rbar)^2;

%% Countinued from 4.4
fun = @(x) besseli(1, x)./besseli(0, x) - rbar;
kappabar = fzero(fun,1);
histogram(data,10,'Normalization','pdf')
hold on
theta = 160:0.1:340;
plot(theta,exp(kappabar*cosd(theta - thetabar))/(360*besseli(0,kappabar)))
ylabel('Probablility Density'); xlabel('Angle (degrees)');
%empirical pdf of paleocurrent data compared to a Von Mises distribution
%with parameters kappabar = 0.26757, mubar = 247.6166

%% The Laws of Large Numbers
x = unidrnd(10,1000,1);
y = [];
for i = 1:1000
    y = [y mean(x(1:i))];
end
plot(y)
title(['The Law of Large Numbers'; 'DiscreteUniform Distribution'])
ylabel('Sample Mean'); xlabel('Sample Size');
hold on
x = 5.5*ones(size(y));
plot(x)
hold off

%% Example 4.5

%% Example 4.6
x = unifrnd(-1, 1, 1000, 1000);
hold on
subplot(2,2,1); histogram(x(1,:))
subplot(2,2,2); histogram(mean(x(1:10,:)))
subplot(2,2,3); histogram(mean(x(1:100,:)))
subplot(2,2,4); histogram(mean(x(1:1000,:))) % CLT applies
hold off

x =  trnd(1, 1000, 1000); %Cauchy distribution is t distribution with 1 dof
hold on
subplot(1,3,1); histogram(x(1,:));
subplot(1,3,2); histogram(mean(x(1:10,:)));
subplot(1,3,3); histogram(mean(x)); %CLT is not increasingly concentrated as N rises, thus 
% the CLT does not apply. This makes sense, because Cauchy distribution mu
% and var is undefined.
hold off

%% Probability Integral Transform
x = unifrnd(0,1,10000,1);
y = raylinv(x, 1);
clf;
h = histogram(y,50,'Normalization','pdf');
hold on
xx = h.BinEdges + h.BinWidth/2;
plot(xx,raylpdf(xx,1),'LineWidth',2)
hold off

%% Example 4.8
x100 = exprnd(1,100,1);
x1000 = exprnd(1,1000,1);
[f1,x1] = ecdf(x100);
[f2,x2] = ecdf(x1000);
plot(x1,f1,x2,f2)
subplot(1,2,1); ecdfhist(f1,x1,20);
subplot(1,2,2); ecdfhist(f2,x2,20);
%empirical cdf is an unbiased estimator of the population cdf

%% Example 4.9
clf;
x = normrnd(0,1,1,5000);
x = [x normrnd(3,1,1,5000)];
[f,xi] = ksdensity(x);
plot(xi,f)
hold on
[f,xi] = ksdensity(x,'bandwidth',0.1);
plot(xi,f,'k--')
[f,xi] = ksdensity(x,'bandwidth',1);
plot(xi,f,'k:')
plot(xi, 0.5*(normpdf(xi,0,1) + normpdf(xi,3,1)),'k-.')
hold off

%% Example 4.10
clf;
x = betarnd(0.5, 0.5, 10000, 1);
[f, xi,bw] = ksdensity(x);
plot(xi,f)
hold on
[f, xi] = ksdensity(x, 'support', [0, 1]);
plot(xi, f, 'k--')
plot(xi, betapdf(xi, 0.5, 0.5), 'k:')
hold off
ylim([0 5]) %necessary to see what is going on
legend('no support','support','arcsine dist')

%% Example 4.12 - stabilized pp plots
x = trnd(1,1,1000);
u = ((1:1000) - 0.5)/1000;
r = 2/pi*asin(sqrt(u));
s = 2/pi*asin(sqrt(normcdf(sort(x))));
subplot(1,2,1); plot(r,s,'k+ ')
axis square
grid on
ylabel 'Transformed Data'; xlabel 'Sine Qunatiles';
title 'Draws from stand. Cauchy dist. using stand. Gaussian';

s = 2/pi*asin(sqrt(tcdf(sort(x), 1)));
subplot(1,2,2); plot(r,s,'k+ ')
axis square
grid on
ylabel 'Transformed Data'; xlabel 'Sine Qunatiles';
title 'Draws from stand. Cauchy dist. using stand. Cauchy';
%straight line

%% Example 4.13 - stabilized pp plots
clf;
x = raylrnd(1,1,1000);
u = ((1:1000) - 0.5)/1000;
r = 2/pi*asin(sqrt(u));
s = 2/pi*asin(sqrt(raylcdf(sort(x))));
subplot(1,3,1); plot(r,s,'k+')
axis square
grid on
ylabel 'Transformed Data'; xlabel 'Sine Qunatiles'; title 'Raleigh Dist.'

pHat = lognfit(x);
s = 2/pi*asin(sqrt(logncdf(sort(x),pHat(1),pHat(2))));
subplot(1,3,2); plot(r,s,'k+')
axis square
grid on
xlabel 'Sine Qunatiles'; title 'Lognormal Dist.'

pHat = gevfit(x); %produces a [k sigma mu]
s = 2/pi*asin(sqrt(gevcdf(sort(x),pHat(1),pHat(2),pHat(3))));
subplot(1,3,3); plot(r,s,'k+')
axis square
grid on
xlabel 'Sine Qunatiles'; title 'Generalized Extreme Value Dist.'

%% Example 4.14 - qq plots
clf;
x = trnd(1,1,10000);
u = ((1:10000) - 0.5)/10000;
subplot(1,2,1); plot(norminv(u), sort(x), 'k+')
axis square
xlabel 'Normal Quantiles'; ylabel 'Sorted Data';
subplot(1,2,2); plot(tinv(u,1),sort(x),'k+')
axis square
xlabel 'Cauchy Quantiles'; ylabel 'Sorted Data';

%% Example 4.15 - qq plots
clf;
n = 1000;
x = raylrnd(1,1,n);
u = ((1:n) - 0.5)/n;
subplot(1,3,3); plot(raylinv(u,1),sort(x),'k+')
axis square
grid on
ylim([0 4]);
xlabel 'Rayleigh Qunatiles';

pHat = gevfit(u);
subplot(1,3,2); plot(gevinv(u,pHat(1),pHat(2),pHat(3)),sort(x),'k+')
axis square
grid on
ylim([0 4]);
xlabel 'GEV Qunatiles';

pHat = lognfit(u);
subplot(1,3,1); plot(logninv(u,pHat(1),pHat(2)),sort(x),'k+')
axis square
grid on
ylim([0 4]);
xlabel 'Log Normal Qunatiles'; ylabel 'Sorted Data';

%% Example 4.16
n = 10000;
x = raylrnd(1, n, 1);
y = normrnd(0, 1, n, 1);
z = x./y;
ksdensity(z,'NumPoints',n) %10000 kernels when n = 10000
xlim([-75 75]) %truncated abscissa
ylabel ' Simulated PDF';

%the ratio of rvs will often result in algebraic tails and infinite
%variance 

%% 4.16 Continued
n = 100;
x = raylrnd(1, n, 1);
y = normrnd(0, 1, n, 1);
z = x.*y;
ksdensity(z,'NumPoints',n) %100 kernels, well behaved
ylabel ' Simulated PDF';

%% Plotting Student's t with different degrees of freedom verse Normal 
clf;
x = -10:0.1:10;

hold on
for nu = 1:5
    plot(x,tpdf(x,nu))
end

pd = makedist('Normal');
plot(x,pdf(pd,x),'k-','LineWidth',2)
ylabel 'Probability Density';
xlim([-4 4])
hold off

legend('1','2','3','4','5','Normal')

%% PDF for correlation coefficent
clf;
r = -1:0.01:1;
N = 10; %degrees of freedom
rho = 0.2; %population correlation coefficent
Pdf = [];

for i = 1:length(r)
    Pdf = [Pdf (N-2)*gamma(N-1)*(1-rho^2)^((N-1)/2)*(1-r(i)^2)^((N-4)/2)*...
        Hypergeometric2f1(1/2, 1/2, N-1/2, (1+rho*r(i))/2)/(sqrt(2*pi)*gamma(N-1/2)*...
        (1-rho*r(i))^(N-3/2))];
end

plot(r,Pdf)

%central moments
fun = @(r, rho, N) r*(N-2)*gamma(N-1)*(1-rho^2)^((N-1)/2)*(1-r^2)^((N-4)/2)*...
        Hypergeometric2f1(1/2, 1/2, N-1/2, (1+rho*r)/2)/(sqrt(2*pi)*gamma(N-1/2)*...
        (1-rho*r)^(N-3/2));

%MATLAB numeric integration
expectedVal = integral(@(r) fun(r, rho, N), -1, 1, 'ArrayValued', true);

%numerical variance calculation
fun1 = @(r, rho, N, xm) (r - xm)^2*(N-2)*gamma(N-1)*(1-rho^2)^((N-1)/2)*...
    (1-r^2)^((N-4)/2)*Hypergeometric2f1(1/2, 1/2, N-1/2, (1+rho*r)/2)/(sqrt(2*pi)*...
    gamma(N-1/2)*(1-rho*r)^(N-3/2));

xm = integral(@(r) fun(r, rho, N), -1, 1, 'ArrayValued', true); %the expected value
Var = integral(@(r) fun1(r, rho, N, xm), -1, 1, 'ArrayValued', true);

%% Order statistics and Example 4.22
clf;
x = 0:0.01:1;
Pdf = zeros(1,length(x));
N = 10;
rlist = [1 3 7 10];

hold on
for k = 1:length(rlist)
    r = rlist(k);
    for i = 1:length(x)
        f = 1/beta(r,N-r+1)*x(i)^(r-1)*(1-x(i))^(N-r);
        Pdf(i) = f;
    end
    plot(x,Pdf)
end
ylim([0 10])
ylabel 'Probability Density';
hold off
%pdf for order statistics of population N whos parent distribution is
%uniform from [0 1]

%expected value
xm = r/(N+1);

%% Previous work continued, pdf of sample median
%median pdf
x = normrnd(0,1,1000,1);
x = sort(x);
Pdf = @(x, N) normcdf(x).^floor(N/2).*(1-normcdf(x)).^(N-floor(N/2)-1)...
    .*normpdf(x)/beta(floor(N/2) + 1, N - floor(N/2));

clf;
hold on
list = [3 11 101];
for i = 1:length(list)
    plot(x,fun(x, list(i)))
end
legend(string(list(1)),string(list(2)),string(list(3)))
xlim([-2 2]); ylim([0 4]); ylabel 'Median PDF'
title 'Sampling dist. of median for stand. norm. parent populations (N) of diff. sizes'

%expected value of the median pdf
xm = integral(@(x) x.*Pdf(x,N), -Inf, Inf); %floating point zero
%variance of the median pdf numerically calculated
xv = integral(@(x) x.^2.*Pdf(x,N), -Inf, Inf);

%% Distribution of the Interquartile Range
clear; clf;
Nlist = [4 12 48 100];

hold on
for i = 1:length(Nlist)
    
    N = Nlist(i);

    fun = @(x, w, N) normcdf(x).^(N/4 - 1).*normpdf(x).*(normcdf(x + w) - normcdf(x)).^...
        (N/2 - 1).*normpdf(x + w).*(1 - normcdf(x + w)).^(N/4);
    f = [];
    s = integral2(@(x, w) fun(x, w, N), -Inf, Inf, 0, Inf);
    
    %analystic integration is intractabe, so normalizing to unity this way
    for w = 0:0.01:5
        f = [f integral(@(x) fun(x, w, N), -Inf, Inf)/s];
    end

    plot(f)
end
ylabel 'IQR PDF';
hold off
%my plot is not the same as Dr. Chave's. I am getting x values that seem
%exactly 100x too large. But I do not know if this is a product of me not
%normalizing correctly or what. Anyways, the shapes are correct

%expected value - choose N and compute
N = Nlist(1);
fun1 = @(x, w, N) w.*normcdf(x).^(N/4 - 1).*normpdf(x).*(normcdf(x + w) - normcdf(x)).^...
    (N/2 - 1).*normpdf(x + w).*(1 - normcdf(x + w)).^(N/4);
e = integral2(@(x, w) fun1(x, w, N), -Inf, Inf, 0, Inf,'AbsTol',1e-12,'RelTol',1e-8);
e = e/integral2(@(x, w) fun(x, w, N), -Inf, Inf, 0, Inf);






