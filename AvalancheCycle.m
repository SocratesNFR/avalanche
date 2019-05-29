% Cycle avalanche analysis

folders = {'MEA2/', 'MEA8/', 'MEA13/', 'MEA14/', ...
    'MEA15/', 'MEA16/', 'MEA17/', 'MEA18/'};
numMEAs = length(folders);
pNLR = cell(numMEAs,1);
pLS = cell(numMEAs,1);
numAvalanches = cell(numMEAs,1);
alphNLR = cell(numMEAs,1);
alphLS = cell(numMEAs,1);
deltaPNLR = cell(numMEAs,1);
deltaPLS = cell(numMEAs,1);

for MEAind = 1:length(folders)
    close all
    % Get list of files for each MEA
    folder = folders{MEAind};
    loc = '/Users/kri/Documents/MATLAB/Avalanche/';
    loc = [loc, folder];
    files = dir(strcat(loc,'*.mat'));
    numFiles = length(files);
    
    % Initialize variables
    pNLR{MEAind} = zeros(numFiles, 1);
    pLS{MEAind} = zeros(numFiles, 1);
    numAvalanches{MEAind} = zeros(numFiles, 1);
    alphNLR{MEAind} = zeros(numFiles, 1);
    alphLS{MEAind} = zeros(numFiles, 1);
    deltaPNLR{MEAind} = zeros(numFiles, 1);
    deltaPLS{MEAind} = zeros(numFiles, 1);
    
    % Avalanche analysis on each file
    for fileind = 1:numFiles
        [pNLR0, pLS0, numAvalanches0, alphNLR0, alphLS0, ...
            deltaPNLR0, deltaPLS0] = avalancheAnalysis(fileind, MEAind);
        pNLR{MEAind}(fileind) = pNLR0;
        pLS{MEAind}(fileind) = pLS0;
        numAvalanches{MEAind}(fileind) = numAvalanches0;
        alphNLR{MEAind}(fileind) = alphNLR0;
        alphLS{MEAind}(fileind) = alphLS0;
        deltaPNLR{MEAind}(fileind) = deltaPNLR0;
        deltaPLS{MEAind}(fileind) = deltaPLS0;
    end
    
end