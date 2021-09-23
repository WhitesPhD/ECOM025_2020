%

clear all; clc; pause(0.01), randn('seed',3212), rand('seed',3212), warning off

Data                    = readtable('StockIndex.xlsx','sheet','WRDS','ReadVariableNames',true);

Returns                 = tick2ret(Data.Index);
SampleSize              = length(Returns);

TestWindowStart         = 10000;
TestWindow              = TestWindowStart : SampleSize;
EstimationWindowSize    = 1000;

pVaR                    = [0.05 0.01];

Historical95            = zeros(length(TestWindow),1);
Historical99            = zeros(length(TestWindow),1);

for t = TestWindow

    i = t - TestWindowStart + 1;
    EstimationWindow = t-EstimationWindowSize:t-1;
    X = Returns(EstimationWindow);
    Historical95(i) = - quantile(X,pVaR(1));
    Historical99(i) = - quantile(X,pVaR(2));

end

figure;
plot(Data.Date(TestWindow),Returns(TestWindow),'linewidth',2),hold on, plot(Data.Date(TestWindow),[-Historical95 -Historical99],'linewidth',2)
ylabel('VaR')
xlabel('Date')
set(gca, 'fontsize', 26)
legend({'Returns','95% Confidence Level','99% Confidence Level'},'Location','Best','AutoUpdate','off')
recessionplot 


% using the riskmetrics

Lambda = 0.94;
Sigma2     = zeros(length(Returns),1);
Sigma2(1)  = Returns(1)^2;

for i = 2 : (TestWindowStart-1)
    Sigma2(i) = (1-Lambda) * Returns(i-1)^2 + Lambda * Sigma2(i-1);
end


Zscore = norminv(pVaR);
EWMA95 = zeros(length(TestWindow),1);
EWMA99 = zeros(length(TestWindow),1);

for t = TestWindow
    k     = t - TestWindowStart + 1;
    Sigma2(t) = (1-Lambda) * Returns(t-1)^2 + Lambda * Sigma2(t-1);
    Sigma = sqrt(Sigma2(t));
    EWMA95(k) = -Zscore(1)*Sigma;
    EWMA99(k) = -Zscore(2)*Sigma;
end

figure;
plot(Data.Date(TestWindow),Returns(TestWindow),'linewidth',2),hold on, plot(Data.Date(TestWindow),[-EWMA95 -EWMA99],'linewidth',2)
ylabel('VaR')
xlabel('Date')
set(gca, 'fontsize', 26)
legend({'Returns','95% Confidence Level','99% Confidence Level'},'Location','Best','AutoUpdate','off')
recessionplot 


% using an asymmetric GARCH

Mdl    = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
EstMdl = estimate(Mdl,Returns(TestWindow));
Sigma2 = infer(EstMdl,Returns(TestWindow));
EGARCH95 = -Zscore(1)*Sigma2.^.5;
EGARCH99 = -Zscore(2)*Sigma2.^.5;


figure;
plot(Data.Date(TestWindow),Returns(TestWindow),'linewidth',2),hold on, plot(Data.Date(TestWindow),[-EGARCH95 -EGARCH99],'linewidth',2)
ylabel('VaR')
xlabel('Date')
set(gca, 'fontsize', 26)
legend({'Returns','95% Confidence Level','99% Confidence Level'},'Location','Best','AutoUpdate','off')
recessionplot 



% zoom-in and compare the VaR Estimates

ZoomInd   = (Data.Date(TestWindow) >= datestr('31-Aug-2007','local')) & (Data.Date(TestWindow) <= datestr('31-Aug-2009','local'));
VaRData   = [-EGARCH95(ZoomInd) -Historical95(ZoomInd) -EWMA95(ZoomInd)  ];
VaRFormat = {'-','--','-.'};

D = Data.Date(TestWindow);
D = D(ZoomInd);

R = Returns(TestWindow);
R = R(ZoomInd);

N = EGARCH95(ZoomInd);
H = Historical95(ZoomInd);
E = EWMA95(ZoomInd);
IndN95    = (R < -N);
IndHS95   = (R < -H);
IndEWMA95 = (R < -E);
figure;
plot(D,R,'LineWidth',1);
hold on
for i = 1 : size(VaRData,2)
    plot(D,VaRData(:,i),VaRFormat{i},'LineWidth',1);
end
ylabel('VaR')
xlabel('Date')
legend({'Returns','EGARCH','Historical','EWMA'},'AutoUpdate','Off')
recessionplot
title('95% VaR violations for different models')
set(gca, 'fontsize', 26)
plot(D(IndN95),-N(IndN95),'o',D(IndHS95),-H(IndHS95),'o',...
   D(IndEWMA95),-E(IndEWMA95),'o','MarkerSize',12,'LineWidth',2)
hold off;




























