% Cycle avalanche analysis
% Run avalanche analysis 

% ########################################################## %
% Include sub-folders containing desired dataset. Folders here are used to
% represent different MEAs in the dataset, and files in the folders are
% recordings from different dates.

% ---- Folders are given as cell arrays of strings
folders = {'', ''};
% ########################################################## %


% Initialize variables
numMEAs = length(folders);
pNLR = cell(numMEAs,1);
pMLE = cell(numMEAs,1);
numAvalanches = cell(numMEAs,1);
fitAvalanches = cell(numMEAs,1);
alphNLR = cell(numMEAs,1);
alphMLE = cell(numMEAs,1);
MFR = cell(numMEAs,1);


for MEAind = 1:length(folders)
    % Get list of files for each MEA
    folder = folders{MEAind};
    
% ########################################################## %    
% Location of dataset folders

    loc = '';
% ########################################################## %

    loc = [loc, folder];
    files = dir(strcat(loc,'*.mat'));
    numFiles = length(files);
    
    % Initialize variables
    pNLR{MEAind} = zeros(numFiles, 1);
    pMLE{MEAind} = zeros(numFiles, 1);
    
    fitAvalanches{MEAind} = zeros(numFiles, 1);
    numAvalanches{MEAind} = zeros(numFiles, 1);
    
    alphNLR{MEAind} = zeros(numFiles, 1);
    alphMLE{MEAind} = zeros(numFiles, 1);
    
    MFR{MEAind} = zeros(numFiles, 1);
    
    % Avalanche analysis on each file
    for fileind = 1:numFiles
        [pNLR0, pMLE0, numAvalanches0, fitAvalanches0, ...
            alphNLR0, alphMLE0, MFR0]...
            = avalancheAnalysis(fileind, MEAind);
        pNLR{MEAind}(fileind) = pNLR0;
        pMLE{MEAind}(fileind) = pMLE0;
        numAvalanches{MEAind}(fileind) = numAvalanches0;
        fitAvalanches{MEAind}(fileind) = fitAvalanches0;
        alphNLR{MEAind}(fileind) = alphNLR0;
        alphMLE{MEAind}(fileind) = alphMLE0;
        MFR{MEAind}(fileind) = MFR0;
    end
    close all
end

savefile = 'results_6sig_last3';
save(savefile, 'alphNLR', 'alphMLE', 'MFR', ...
    'numAvalanches', 'fitAvalanches', 'pNLR', 'pMLE');