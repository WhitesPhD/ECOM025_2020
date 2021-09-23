%% Week 3 Financial Econometrics
% Daniele Bianchi
% d.bianchi@qmul.ac.uk

%--------------------------------------------------------------------------
% Housekeeping
%--------------------------------------------------------------------------
clear all; clc; pause(0.01), randn('seed',3212), rand('seed',3212), warning off

%--------------------------------------------------------------------------
% Load data from yahoo finance
%--------------------------------------------------------------------------

Data            = readtable('StockIndex.xlsx','Sheet','Multi','ReadVariableNames',true);

%Plot some properties of the returns on stocks and bonds

figure
plot(Data.Date,Data.StockMkt,'linewidth',2)
set(gca, 'fontsize', 26)
legend({'Stock Market Returns'},'AutoUpdate','off')
recessionplot 


figure
plot(Data.Date,Data.BondReturns,'linewidth',2)
set(gca, 'fontsize', 26)
legend({'Bond Market Returns'},'AutoUpdate','off')
recessionplot 

nbins = 10;
figure
subplot(1,2,1)
histogram(Data.StockMkt,nbins),title('Stock Market')
set(gca, 'fontsize', 26)
subplot(1,2,2)
histogram(Data.BondReturns,nbins),title('Bond Market')
set(gca, 'fontsize', 26)

figure
subplot(2,1,1)
qqplot(Data.StockMkt),title('Stock Market')
set(gca, 'fontsize', 26)
subplot(2,1,2)
qqplot(Data.BondReturns),title('Bond Market')
set(gca, 'fontsize', 26)

% calculate and display the correlation between stocks and bonds; are bonds
% good for hedging or for diversification?

y = [Data.StockMkt Data.BondReturns];
rho  = corr(y);
display(['The correlation between stocks and bond is: ', num2str(rho(1,2))])

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
plot(Data.Date,EWMArho,'linewidth',2),ylim([-0.7 0.7])
set(gca, 'fontsize', 26)
legend({'Correlation between Stock and Bond Market (EWMA)'},'AutoUpdate','off')
recessionplot 
 

[p, lik, Ht]  = dcc_mvgarch(y,1,1,1,1);
Ht = reshape(Ht,4,T)';
DCCrho = Ht(:,3) ./ sqrt(Ht(:,1) .* Ht(:,4)); %% DCCrho is a vector of correlations

figure
plot(Data.Date,DCCrho,'linewidth',2),ylim([-0.3 0.45])
set(gca, 'fontsize', 26)
legend({'Correlation between Stock and Bond Market (EWMA)'},'AutoUpdate','off')
recessionplot 

% superimpose the two graphs

figure
plot(Data.Date,EWMArho,'linewidth',2),ylim([-0.5 0.65])
hold on
plot(Data.Date,DCCrho,'linewidth',2)
set(gca, 'fontsize', 26)
legend({'EWMA','DCC'},'AutoUpdate','off')
recessionplot 


Data.sqret = Data.returns.^2;

ma_1 = tsmovavg(Data.sqret','s',90);
ma_2 = tsmovavg(Data.sqret','s',180);
ma_3 = tsmovavg(Data.sqret','s',720);

figure
plot(Data.Date,ma_1)
hold on
plot(Data.Date,ma_2)
hold on
plot(Data.Date,ma_3)
set(gca, 'fontsize', 26)
legend({'m1','ma2','ma3'},'AutoUpdate','off')
recessionplot 

% qqplot
figure
qqplot(Data.returns)
set(gca, 'fontsize', 26)
title('qq-plot standard normal')

% estimate a simple GARCH(1,1)

Mdl     = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,Data.returns);
v       = infer(EstMdl,Data.returns);

% plot the returns

figure
plot(Data.Date,Data.returns,'linewidth',1),ylim([-25 12])
set(gca, 'fontsize', 26)
legend({'Daily returns'},'AutoUpdate','off')
recessionplot 

% plot the conditional volatility 

figure
plot(Data.Date,sqrt(v),'linewidth',1)
set(gca, 'fontsize', 26)
legend({'Estimated volatility'},'AutoUpdate','off')
recessionplot 

% plot returns and conditional volatility

lb = - 2*sqrt(v);
ub = 2*sqrt(v);

figure
plot(Data.Date,Data.returns,'linewidth',1),title('Returns and conditional volatility')
hold on
plot(Data.Date,[lb ub],'color','red','linewidth',1)
set(gca, 'fontsize', 26)
recessionplot 

% calculate the standardized residuals 

z = Data.returns./sqrt(v);

figure
subplot(1,2,1)
qqplot(Data.returns)
set(gca, 'fontsize', 26)
title('qq-plot returns')
subplot(1,2,2)
qqplot(z)
set(gca, 'fontsize', 26)
title('qq-plot standardized returns')

% autocorrelation of returns vs standardized returns 

nLags = 60;
[acf, ~, bounds]   = autocorr(Data.returns, nLags) ; 
[acfs, ~, boundss] = autocorr(z, nLags) ; 

figure
subplot(2,1,1)
plot(acf(2 : end), 'o', 'MarkerFaceColor', 'blue','linewidth',2),ylim([-0.1 0.1]),title('ACF returns')
hold on
plot(ones(nLags, 1)*bounds', '-.', 'color', 'red','linewidth',2)
set(gca, 'fontsize', 26)
subplot(2,1,2)
plot(acfs(2 : end), 'o', 'MarkerFaceColor', 'blue','linewidth',2),ylim([-0.1 0.1]),title('ACF standardized returns')
hold on
plot(ones(nLags, 1)*boundss', '-.', 'color', 'red','linewidth',2)
set(gca, 'fontsize', 26)

% residuals analysis
% notice the null hypothesis that the data comes from a normal distribution with an unknown mean and variance
% The result h is 1 if the test rejects the null hypothesis at the 5% significance level, and 0 otherwise.

[h,p,jbstat,critval] = jbtest(z);

% estimate a GARCH(1,1) with Student-t

Mdl_t     = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN,'Distribution','t');
[EstMdl_t,EstParamCov_t,logL_t,info_t] = estimate(Mdl_t,Data.returns);

v_t       = infer(EstMdl_t,Data.returns);

figure
plot(Data.Date,sqrt(v),'linewidth',1)
set(gca, 'fontsize', 26)
hold on
plot(Data.Date,sqrt(v_t),'linewidth',1)
set(gca, 'fontsize', 26)
legend({'Gaussian','Student-t'},'AutoUpdate','off')
recessionplot 


z = Data.returns./sqrt(v);

[h,p,jbstat,critval] = jbtest(z);


% Likelihood ratio between a Normal GARCH and a Student-t GARCH
% returns a logical value (h) and the pvalue of the rejection decision from conducting a likelihood ratio test of model specification.
% the null hypothesis of the LR is that the models are the same
% here the model with Student-t is the unrestricted one

q = 1; %number of restrictions
[h,pValue,stat,cValue] = lratiotest(logL_t,logL,q);





































