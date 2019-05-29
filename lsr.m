% --- Least squares regression for log-log plot of avalanche size
% Fit linear regression to log(frequency) vs. log(avalanche size) plot
% Exclude first bin and final bins with total probability <1‰ that of the
% first bin, following Pasquale et al. (2008) and Massobrio et al. (2015)
% Excluded only final bin: 60 electrodes

function [pdfLS, alphLS, genZetaLS] = lsr(xmin, xmax, pdfEMP)

x = xmin:xmax;
y = pdfEMP;
X = log(x).';
Y = log(y).';
X = X(y>0);
Y = Y(y>0);
[LSfit, gof] = fit(X, Y, 'poly1');

% Least-squares fitting PDF
coeffsLS = coeffvalues(LSfit);
alphLS = -coeffsLS(1);
genZetaLS = sum(((0:1000000) + xmin).^(-alphLS));
pdfLS = x.^(-alphLS)/genZetaLS;