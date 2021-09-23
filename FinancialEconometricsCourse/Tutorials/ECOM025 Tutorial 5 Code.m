%% Week 5 Financial Econometrics
% Daniele Bianchi
% d.bianchi@qmul.ac.uk

%--------------------------------------------------------------------------
% Housekeeping
%--------------------------------------------------------------------------
clear all; clc; pause(0.01), randn('seed',3212), rand('seed',3212), warning off

%--------------------------------------------------------------------------
% Load data from an excel spreadsheet
%--------------------------------------------------------------------------

Data            = readtable('DataFX.xls','ReadVariableNames',true);

%Calculate and plot the returns

FX      =   [Data.USDEUR Data.USDGBP Data.USDAUS];
Returns =   diff(log(FX))*100;


figure
subplot(3,1,1)
plot(Data.Date(2:end),Returns(:,1),'linewidth',2)
set(gca, 'fontsize', 26)
legend({'Returns on the USD/EUR'},'AutoUpdate','off')
recessionplot 

subplot(3,1,2)
plot(Data.Date(2:end),Returns(:,2),'linewidth',2)
set(gca, 'fontsize', 26)
legend({'Returns on the USD/GBP'},'AutoUpdate','off')
recessionplot 

subplot(3,1,3)
plot(Data.Date(2:end),Returns(:,3),'linewidth',2)
set(gca, 'fontsize', 26)
legend({'Returns on the USD/AUS'},'AutoUpdate','off')
recessionplot 


% statistical properties of the returns on the USD/GBP

nbins = 10;
figure
subplot(2,1,1)
histogram(Returns(:,2),nbins),title('Returns on the USD/EUR')
set(gca, 'fontsize', 26)

subplot(2,1,2)
qqplot(Returns(:,2)),title('Returns on the USD/EUR')
set(gca, 'fontsize', 26)


% estimate a simple GARCH(1,1)

Mdl     = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,Returns(:,2));
v       = infer(EstMdl,Returns(:,2));

% plot the volatility

figure
plot(Data.Date(2:end),sqrt(v),'linewidth',1)
set(gca, 'fontsize', 26)
legend({'GARCH(1,1) Volatility'},'AutoUpdate','off')
recessionplot 



% Calculate the unconditional correlation between USD/GBP and USD/EUR

y    = Returns(:,1:2);
rho  = corr(y);
display(['The correlation between USD/GBP and USD/EUR: ', num2str(rho(1,2))])

% dynamic correlation from an EWMA

T           = length(y);
EWMA        = nan(T,3); % create a matrix to hold covariance matrix for each t
lambda      = 0.94;
S           = cov(y); % initial (t=1) covar matrix 
EWMA(1,:)   = S([1,4,2]); % extract var and covar

for i = 2:T % loop though the sample
S = lambda*S+(1-lambda)* y(i-1,:)'*y(i-1,:); EWMA(i,:) = S([1,4,2]); % convert matrix to vector
end
EWMArho = EWMA(:,3)./sqrt(EWMA(:,1).*EWMA(:,2)); % calculate correlations


figure
plot(Data.Date(2:end),EWMArho,'linewidth',2)
set(gca, 'fontsize', 26)
legend({'Correlation between USD/GBP and USD/EUR (EWMA)'},'AutoUpdate','off')
recessionplot 
 



































