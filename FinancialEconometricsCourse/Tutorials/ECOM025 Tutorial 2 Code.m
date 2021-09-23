%% Week 1 Financial Econometrics
% Daniele Bianchi
% d.bianchi@qmul.ac.uk

%--------------------------------------------------------------------------
% Housekeeping
%--------------------------------------------------------------------------
clear all; clc; pause(0.01), randn('seed',3212), rand('seed',3212), warning off

%--------------------------------------------------------------------------
% Load data from yahoo finance
%--------------------------------------------------------------------------

Data            = readtable('StockIndex.xlsx','Sheet','WRDS','ReadVariableNames',true);

returns         = diff(log(Data.Index))*100;
Data.returns    = [0; returns];

% plot the index and the returns on the S&P500 index

figure
subplot(2,1,1),
plot(Data.Date,Data.Index)
set(gca, 'fontsize', 26)
recessionplot
subplot(2,1,2),
plot(Data.Date,Data.returns)
set(gca, 'fontsize', 26)
recessionplot

% define year variables

Data.year = year(Data.Date);
G         = findgroups(Data.year);
vol       = splitapply(@std,Data.returns,G);
date      = splitapply(@mean,Data.year,G);

figure
bar(date, vol*sqrt(250))
set(gca, 'fontsize', 26)


% autocorrelation 

figure
subplot(1,2,1)
autocorr(Data.returns, 360)
subplot(1,2,2)
autocorr(Data.returns.^2, 360)

% testing for autocorrelations

lags        = [20 180 360];
pvalue_r    = [];
pvalue_s    = [];
stat_r      = [];
stat_s      = [];

for i=1:3

[h,pValue,stat]         = lbqtest(Data.returns,'lags',lags(i));
stat_r   = [stat_r;stat];
pvalue_r = [pvalue_r;pValue];

[h,pValue,stat]   = lbqtest(Data.returns.^2,'lags',lags(i));
stat_s   = [stat_s;stat];
pvalue_s = [pvalue_s;pValue];

end


% empirical density vs normal

figure
histfit(Data.returns,100,'normal')


% empirical vs theoretical CDF

rng('default')  % For reproducibility

cdfplot(Data.returns)
hold on
x = linspace(min(Data.returns),max(Data.returns));
plot(x,normcdf(x,mean(Data.returns),std(Data.returns)))
legend('Empirical CDF','Theoretical CDF','Location','best')
hold off

% qqplot
figure
qqplot(Data.returns)

pd = makedist('tLocationScale','mu',mean(Data.returns),'sigma',std(Data.returns),'nu',3);

figure
qqplot(Data.returns,pd)



