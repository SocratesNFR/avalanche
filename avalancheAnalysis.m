%% Neuronal avalanches
% Reads data from event detection files, computes the average interspike
% interval and determines sizes of avalanches based on the definition given 
% by Massobrio et al. (2015). Computes p-value based on synthetic data
% generated from hypothesized power-law PDF using the method by Clauset et
% al. (2009).

function [pNLR, pLS, totAvalanches, alphNLR, alphLS, deltaPNLR, deltaPLS]...
    = avalancheAnalysis(fileind, MEAind)

%% Clean up
%clc; clear;

%% Read event file
% rate = 10e3;
%fileind = 14;
%MEAind = 1;

folders = {'MEA2/', 'MEA8/', 'MEA13/', 'MEA14/', ...
    'MEA15/', 'MEA16/', 'MEA17/', 'MEA18/'};
folder = folders{MEAind};

loc = '/Users/kri/Documents/MATLAB/Avalanche/5sig/';
loc = [loc, folder];
files = dir(strcat(loc,'*.mat'));
filename = files(fileind).name;
load(strcat(loc,filename));
evs = load(strcat(loc,filename), 'events');
evs = evs.events;

%% Noisy electrodes
excl_list = cell(8, 18);
excl_list{2, 3} = 42;       % MEA 8, 05.09, electrode 64
excl_list{3, 3} = 42;       % MEA 13, 05.09, electrode 64
excl_list{3, 17} = 42;      % MEA 13, 02.10, electrode 64
excl_list{5, 5} = 42;       % MEA 15, 11.09, electrode 64
excl_list{6, 3} = 42;       % MEA 16, 05.09, electrode 64
excl_list{6, 5} = [42 24];  % MEA 16, 11.09, electrodes 64, 42
excl_list{6, 12} = 10;      % MEA 16, 26.09, electrode 24
excl_list{6, 14} = 42;      % MEA 16, 29.09, electrode 64
excl_list{7, 18} = 42;      % MEA 17, 03.10, electrode 64
% for m = 1:8
%     for n = 1:18
%         excl_list{m,n} = [excl_list{m,n}, 42];
%     end
% end
excl = excl_list{MEAind, fileind};

%% Plot titles
MEAnum = {'MEA2', 'MEA8', 'MEA13', 'MEA14', ...
    'MEA15', 'MEA16', 'MEA17', 'MEA18'};
recDate = {filename(end-8:end-4)};
plotTitle = strcat(MEAnum{MEAind}, {', Date: '}, recDate);

%% Avalanche detection
[totEvents, avgISI, avalancheID, avalancheSize, numAvalanches] = ...
    avalancheDetection(evs, time, rate, excl);

%% Get histogram counts for avalanche sizes
sizeBins = 1:60;
edges = 0.5:1:60.5;
[sizeFreqs, edges] = histcounts(avalancheSize, edges);
totAvalanches = sum(sizeFreqs);

% Empirical PDF
xmin = 2;
xmax = 59;
x = xmin:xmax;
nS = sum(sizeFreqs(x));
pdfEMP = sizeFreqs(x)/nS;

%% Fitting power law for estimate of coefficients
[pdfLS, alphLS, genZetaLS] = lsr(xmin, xmax, pdfEMP);
% Delta p defined according to eq. (11) of Tetzlaff et al. (2010) as the 
% mean of the differences between empirical and fitted PDFs
deltaPLS = mean(pdfEMP - pdfLS);

%% Nonlinear regression based on initial estimates from LS
[pdfNLR, alphNLR] = nlr(alphLS, 1/genZetaLS, x, pdfEMP);
deltaPNLR = mean(pdfEMP - pdfNLR);

%% CDFs and Kolmogorov-Smirnov distances
[cdfNLR, cdfEMP, dKSempNLR] = cdfKSdist(pdfNLR, pdfEMP);
[cdfLS, cdfEMP, dKSempLS] = cdfKSdist(pdfLS, pdfEMP);

%% Goodness of fit
pNLR = gofClauset(cdfNLR, dKSempNLR, nS, xmin, xmax);
pLS = gofClauset(cdfLS, dKSempLS, nS, xmin, xmax);

%% Plot results with GoF
figure
loglog(1:60, sizeFreqs, x, pdfNLR*nS);
xlabel('Avalanche size (no. of electrodes)')
ylabel('Frequency')
alphaLeg = string(round(alphNLR,3));
pLeg = string(pNLR);
nLeg = string(totAvalanches);
fitLeg = strcat({'Fit: alpha = '}, alphaLeg, {', p = '}, pLeg);
expLeg = strcat({'Experimental data, N = '}, nLeg);
legend(expLeg, fitLeg)
title(plotTitle)
set(gcf, 'color', 'white')

savefig(strcat(filename(1:end-4), '-NLR'))
saveas(gcf, strcat(filename(1:end-4), '-NLR'), 'png')

figure
loglog(1:60, sizeFreqs, x, pdfLS*nS);
xlabel('Avalanche size (no. of electrodes)')
ylabel('Frequency')
alphaLeg = string(round(alphLS,3));
pLeg = string(pLS);
nLeg = string(totAvalanches);
fitLeg = strcat({'Fit: alpha = '}, alphaLeg, {', p = '}, pLeg);
expLeg = strcat({'Experimental data, N = '}, nLeg);
legend(expLeg, fitLeg)
title(plotTitle)
set(gcf, 'color', 'white')

savefig(strcat(filename(1:end-4), '-LS'))
saveas(gcf, strcat(filename(1:end-4), '-LS'), 'png')