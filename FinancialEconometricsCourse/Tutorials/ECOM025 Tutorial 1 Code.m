% Week 1 Financial Econometrics
% Daniele Bianchi
% d.bianchi@qmul.ac.uk

%--------------------------------------------------------------------------
% Housekeeping
%--------------------------------------------------------------------------
clear all; clc; pause(0.01), randn('seed',3212), rand('seed',3212), warning off

%--------------------------------------------------------------------------
% Load data from the excel file
%--------------------------------------------------------------------------
Data    = readtable('StockIndex.xlsx','Sheet','Multi','ReadVariableNames',true);

% Calculating simple moments

y       = [Data.StockMkt Data.BondReturns];

mean(y)         % average returns
var(y)          % variance of the returns
std(y)          % volatility (square root of the variance)
skewness(y)     % skewness of the returns
kurtosis(y)     % kurtosis of the returns
cov(y)          % covariance between stock and bonds
corr(y)         % correlation between stock and bonds


% plot the time series of monthly stock excess returns
figure
plot(Data.Date,Data.StockMkt,'linewidth',2)
set(gca, 'fontsize', 26)
legend({'Stock Market Returns'},'AutoUpdate','off')
recessionplot


% plot the time series of monthly bond excess returns
figure
plot(Data.Date,Data.BondReturns,'linewidth',2)
set(gca, 'fontsize', 26)
legend({'Bond Returns'},'AutoUpdate','off')
recessionplot

%plot the probability distrubution function of stock market returns
nbins = 20; % define the number of bins in the histrogram

figure
hist(Data.StockMkt,nbins)
set(gca, 'fontsize', 26)

%plot the probability distrubution function of bond returns
nbins = 20; % define the number of bins in the histrogram

figure
hist(Data.BondReturns,nbins)
set(gca, 'fontsize', 26)

% simulate and plot distribution functions
mu      = mean(Data.StockMkt);
sigma   = std(Data.StockMkt);
T       = length(Data.StockMkt);
nu      = 5;

x_norm  = normrnd(mu,sigma,T,1);
x_t     = mu + sigma*trnd(nu,T,1);

nbins = 10; % define the number of bins in the histrogram

figure
subplot(1,2,1),hist(Data.StockMkt,nbins),title('Empirical distribution')
set(gca, 'fontsize', 26)
subplot(1,2,2),hist(x_norm,nbins),title('Simulated distribution (normal)')
set(gca, 'fontsize', 26)

figure
subplot(1,2,1),hist(Data.StockMkt,nbins),title('Empirical distribution')
set(gca, 'fontsize', 26)
subplot(1,2,2),hist(x_t,nbins),title('Simulated distribution (Student-t)')
set(gca, 'fontsize', 26)




