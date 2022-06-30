%% Example 10.10
% An empirical orthogonal function (eof) analysis will be applied to the
% eight zonal bands of climate change data.

%% Importing Data
clear; clf; close all;
z = importdata("ZonAnnTs+dSST.dat");
x = z(:, 8:15);
varx = var(x);

% The variance range covers an orfer of magnitude, with the largest
% variance being near the poles.
% Therefore, the covariance eof's will likely be dominated by those
% variables.

%% Empircal Orthogonal Functions
x1 = x - repmat(mean(x), length(x), 1);

% Singular-Value Decomposition
[u, s, v] = svd(x1);
percvar = diag(s.^2)'/sum(diag(s.^2));

% The first two eof's account for ~87% of the total variance.

%% Scree Plot for EOF
figure(1);
plot(diag(s.^2), 'b-', 'LineWidth', 1.75)

ylim([0 2*10^6])
ylabel 'Squared Singular Value'
xlabel 'EOF'
title 'Scree Plot for Zonally Averaged Temperature Data'

%% Data Analysis
% Examine the eof's:
disp(v)

figure(2);
eofa = u*s; % eof amplitude
t = z(:, 1);
plot(t, eofa(:, 1), 'c-', 'LineWidth', 1.5)
hold on
plot(t, -eofa(:, 2), 'm-', 'LineWidth', 1.5)
% reversing the sign to match with data - eofa's are insenstive to sign
% changes, so examination of actual trends is necessary
hold off

xlim([1880 2020]); ylim([-200 300]);
xlabel 'Year'; ylabel 'EOF Amplitude';
title 'First Two EOF Amplitudes for Unscaled Temperature Data'
legend 'First EOF' 'Second EOF' Location nw

%% Varimax Rotation to First Four Eigenvectors
% Varimax rotation maximizes the sum of squares of eigenvectors. This is
% accomplished by simplifying the columns of the truncated matrix of right
% singular vectors, V.
v1 = rotatefactors(v(:, 1:4), 'Method', 'varimax');

% The varimax rotation results in the north and south poles dominating the
% first two eigenvectors.

figure(3);
plot(t, x*v1(:,1), 'LineWidth', 1.1)
hold on
plot(t, -x*v1(:,2), 'LineWidth', 1.1)
% sign change appied
plot(t, -x*v1(:,3), 'LineWidth', 1.1)
% sign change applied
plot(t, x*v1(:,4), 'LineWidth', 1.1)
hold off

legend 'First Eigenvector' 'Second Eigenvector' 'Third Eigenvector' 'Fourth Eigenvector' Location se
xlabel 'Year'; ylabel 'EOF Amplitude';
title 'First Four EOF Amplitudes After Varimax Rotation'

%% Repeat For Correlation EOF's
% The exercise can be repeated using the correlation eof's rather than the
% covariance eof's

b = diag(1./std(x));
[u, s, v] = svd(x1*b);
percvar = diag(s.^2)'/sum(diag(s.^2));

%% Scree Plot for Data with Columns Normalized by their Standard Deviation
figure(4);
plot(diag(s.^2), 'b-', 'LineWidth', 1.75)

ylim([0 800])
ylabel 'Squared Singular Value'
xlabel 'EOF'
title(["Scree Plot for Zonally Averaged Temperature Data";...
    "with Columns Normalized by Standard Deviations"])

% The scree plot suggests that only two eof's should be retained.

%% Examining the Data
disp(v)

% Plotting EOF Amplitudes
figure(5)
eofa = u*s;
plot(t, eofa(:, 1), 'b-', 'LineWidth', 1.5)
hold on 
plot(t, -eofa(:, 2), 'r-', 'LineWidth', 1.5)
plot(t, eofa(:, 3), 'y-', 'LineWidth', 1.5)
hold off

xlabel 'Year'; ylabel 'EOF Amplitude';
title 'First Three EOF Amplitudes for Normalized Temperature Data'
legend 'First EOF' 'Second EOF' Location nw

% Note the much smaller ordinate for the amplitude measurement in the
% normalized plots compared to the non0normalized plots.

%% Apply Varimax Rotation
v1 = rotatefactors(v(:, 1:3), "Method", "varimax");

% Plotting Varimax Rotated EOF Amplitudes
figure(6)
plot(t, x*b*v1(:, 1), 'b-', 'LineWidth', 1.5)
hold on
plot(t, -x*b*v1(:, 2), 'r-', 'LineWidth', 1.5)
% reversing the sign to match the data
plot(t, x*b*v1(:, 3), 'y-', 'LineWidth', 1.5)
hold off

ylim([-4 6])
xlabel 'Year'; ylabel 'EOF Amplitude';
title 'First Three Varimax Rotated EOF Amplitudes for Normalized Data'
legend 'First Amplitude' 'Second Amplitude' 'Third Amplitude' Location nw

%% Apply Regression Preweighting
% A third version of preweighting based on regressing a given column
% against the remaining ones will now be employed.
p = 8;
res = zeros(1, p);
for i = 1:p
    y = x1(:, i);
    if i == 1
        x2 = x1(:, 2:p);
    elseif i < p
        x2 = [x1(:, 1:(i - 1)) x1(:, i + 1:p)];
    else
        x2 = x1(:, 1:(p - 1));
    end
    b1 = x2\y;
    res(i) = std(y - x2*b1);
end

b = diag(1./res);
varx = var(x1*b);

% This form of weighting emphasizes the equator zones relative to the polar
% ones and hence had the opposite sense to the covariance approach.

[u, s, v] = svd(x1*b);
eofa = u*s;
percvar = diag(s.^2/sum(diag(s.^2)));

%% Scree Plot for Regression Normalized Data
figure(7);
plot(diag(s.^2), 'k-', 'LineWidth', 1.5)

ylim([0 8000])
xlabel 'EOF'; ylabel 'Squared Singular Values'
title 'Scree Plot for Data Normalized by Noise Variance'

%% Examine the Data
disp(v)

% Plot First Three EOF Amplitudes for Data Normalized by Noise Variance
figure(8);
plot(t, x1*b*v(:,1), 'r-', 'LineWidth', 1.5)
hold on
plot(t, x1*b*v(:,2), 'g-', 'LineWidth', 1.5)
plot(t, -x1*b*v(:,3), 'k-', 'LineWidth', 1.5)
% reversing the sign
hold off

xlabel 'Year'; ylabel 'EOF Amplitude';
ylim([-20 30])
title("First Three EOF's for the Data Normalized by Noise Variance")
legend 'First EOF' 'Second EOF' 'Third EOF' Location se

%% Apply Varimax Rotation
v1 = rotatefactors(v(:, 1:3));

figure(9);
plot(t, x1*b*v1(:,1), 'b-', 'LineWidth', 1.5)
hold on
plot(t, x1*b*v1(:,2), 'c-', 'LineWidth', 1.5)
plot(t, -x1*b*v1(:,3), 'm-', 'LineWidth', 1.5)
% reversing the sign
hold off

xlabel 'Year'; ylabel 'EOF Amplitude';
ylim([-20 20])
title("First Three Varimax Rotated EOF's for Data Normalized by Noise Variance")
legend 'First EOF' 'Second EOF' 'Third EOF' Location se