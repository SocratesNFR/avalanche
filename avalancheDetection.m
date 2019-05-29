% --- Performs avalanche detection and computes size distribution from
% event list

function [totEvents, avgISI, avalancheID, avalancheSize, numAvalanches] = ...
    avalancheDetection(events, time, rate, excl)

%% Calculate average interspike interval (network level)
% Event times in s, isi in ms
evTimes = [];
for el = 1:length(events)
    if ~ismember(el,excl) && events{el}(1,3) ~= 0
        evTimes = [evTimes; events{el}(:,2)];
    end
end
totEvents = length(evTimes);
unique_ev_times = unique(evTimes);
unique_ev_times = sort(unique_ev_times);
ISI = diff(1e3*unique_ev_times);
avgISI = mean(ISI <= 100);

%% Time binning with bin width = avgISI

binWidth = round(avgISI*1e-3*rate);
if binWidth == 0
    binWidth = 1;
end
binnedSize = floor(length(time)/binWidth);
spikeTrains = zeros(length(time), length(events));
spikesBinned = zeros(binnedSize, length(events));
for el = 1:length(events)
    % Generate spike trains where 0 = no event at that time, 1 = event at
    % that time
    if ~ismember(el,excl) && events{el}(1,3) ~= 0
        inds = events{el}(:,3);
        spikeTrains(inds,el) = 1;
    end
    % Binned spike counts
    for i = 1:binnedSize
        iStart = (i-1)*binWidth + 1;
        iEnd = i*binWidth;
        spikesBinned(i,el) = sum(spikeTrains(iStart:iEnd,el));
    end
end

%% Avalanche detection
% Create array where 0 means no avalanche and every time bin included in an
% avalanche is labeled with an avalanche ID number
avalancheID = zeros(binnedSize,1);
k = 0;
isAvalanche = 0;
popBinned = sum(spikesBinned,2);
for i = 1:binnedSize
    if popBinned(i) > 0 && isAvalanche == 0
        k = k+1;
        avalancheID(i) = k;
        isAvalanche = 1;
    elseif popBinned(i) > 0 && isAvalanche == 1
        avalancheID(i) = k;
    elseif popBinned(i) == 0 && isAvalanche == 1
        isAvalanche = 0;
    end
end

% Compute size and duration of each avalanche
avalancheSize = zeros(k,1);
avalancheDur = zeros(k,1);
for j = 1:k
    avalancheInds = find(avalancheID == j);
    avalancheDur(j) = length(avalancheInds);   % num bins        %*binWidth/rate*1e3;   % ms
    for el = 1:length(events)
        if sum(spikesBinned(avalancheInds,el)) > 0
            avalancheSize(j) = avalancheSize(j)+1;
        end
    end
end
numAvalanches = length(avalancheSize);