%% Neuronal avalanches
% Reads data from event detection files, computes the average interspike
% interval and determines sizes of avalanches based on the definition given 
% by Massobrio et al. (2015). Computes p-value based on synthetic data
% generated from hypothesized power-law PDF using the method by Clauset et
% al. (2009).

function [pNLR, pMLE, totAvalanches, fitAvalanches, alphNLR, alphMLE, MFR]...
    = avalancheAnalysis(fileind, MEAind)

%% Read event file

% ########################################################## %
% Info about desired dataset
% folders: The sub-folders containing the data to be analyzed
% loc: The path to the sub-folders indicated by folders
% saveloc: The location where the figures will be saved

folders = {'', ''};
folder = folders{MEAind};
loc = '';
saveloc = '';
% ########################################################## %

loc = [loc, folder];
files = dir(strcat(loc,'*.mat'));
filename = files(fileind).name;
load(strcat(loc,filename));
evs = load(strcat(loc,filename), 'events');
evs = evs.events;

%% Noisy electrodes

% ########################################################## %
% List any noisy electrodes to be excluded from the analysis. Noisy
% electrodes can be defined per MEA per recording.

numRecs = length(files);
numMEAs = length(folders);
excl_list = cell(numMEAs, numRecs);
% ########################################################## %

excl = excl_list{MEAind, fileind};

%% Plot titles

% ########################################################## %
% Define plot titles as desired. The titles as suggested here are the MEA
% identifier and the recording date from the filename (from MultiChannel
% Systems format for file naming).

% ---- Primary cortical June 2019
MEAnum = cell(length(folders),1);
for i = 1:length(folders)
    MEAnum{i} = folders{i}(1:end-1);
end
recDate = {filename(1:10)};
plotTitle = strcat(MEAnum{MEAind}, {', Date: '}, recDate);
% ########################################################## %


%% MFR
MFR = mean(FR(FR>0));

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
if sizeFreqs(xmax) == 0
    notEmpty = find(sizeFreqs ~=0);
    xmax = notEmpty(end);
end
x = xmin:xmax;
nS = sum(sizeFreqs(x));
pdfEMP = sizeFreqs(x)/nS;
fitAvalanches = sum(sizeFreqs(x));

%% Nonlinear regression estimate for power law fit
% Least-squares fitting of power law for estimate of coefficients
[pdfLS, alphLS, coeffLS] = lsr(xmin, xmax, pdfEMP);

% Nonlinear regression based on initial estimates from LS
[pdfNLR, alphNLR] = nlr(alphLS, coeffLS, x, pdfEMP);

%% Maximum likelihood estimate for power law fit
powerlawPDF = @(data, alpha) (alpha - 1)*data.^(-alpha);
alphMLE = mle(avalancheSize(avalancheSize>=xmin & avalancheSize<=xmax),...
    'pdf', powerlawPDF, 'start', 1.01);
pdfMLE = x.^(-alphMLE);
pdfMLE = pdfMLE./sum(pdfMLE);

%% CDFs and Kolmogorov-Smirnov distances
[cdfNLR, cdfEMP, dKSempNLR] = cdfKSdist(pdfNLR, pdfEMP);
[cdfMLE, cdfEMP, dKSempMLE] = cdfKSdist(pdfMLE, pdfEMP);

%% Goodness of fit
pNLR = gofClauset(cdfNLR, dKSempNLR, nS, xmin, xmax);
pMLE = gofClauset(cdfMLE, dKSempMLE, nS, xmin, xmax);

%% Plot results with GoF and save

fNLR = figure;
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
set(fNLR, 'color', 'white')
yl = ylim;
ylim([1, yl(2)])
xlim([1, 100])

savefig(fNLR, strcat(saveloc, filename(1:end-4), '-NLR'))
saveas(fNLR, strcat(saveloc, filename(1:end-4), '-NLR'), 'png')


fMLE = figure;
loglog(1:60, sizeFreqs, x, pdfMLE*nS);
xlabel('Avalanche size (no. of electrodes)')
ylabel('Frequency')
alphaLeg = string(round(alphMLE,3));
pLeg = string(pMLE);
nLeg = string(totAvalanches);
fitLeg = strcat({'Fit: alpha = '}, alphaLeg, {', p = '}, pLeg);
expLeg = strcat({'Experimental data, N = '}, nLeg);
legend(expLeg, fitLeg)
title(plotTitle)
set(fMLE, 'color', 'white')
yl = ylim;
ylim([1, yl(2)])
xlim([1, 100])

savefig(fMLE, strcat(saveloc, filename(1:end-4), '-MLE'))
saveas(fMLE, strcat(saveloc, filename(1:end-4), '-MLE'), 'png')