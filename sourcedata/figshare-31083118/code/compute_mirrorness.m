clear binEdges binIdx_all Nbins electrodesInBin propMat binCenters binLabels
clear minEdge maxEdge side_all

% AGGREGATED “PROP MAT” ACROSS MULTIPLE SUBJECTS (DATA‐DRIVEN BINS)
% Loops through subject folders in the current directory,
% loads each subject’s metadata.mat, condensed_ecog_info_only.mat, stimset.mat,
% and accumulates, for each data‐driven olf_coor bin, across all subjects:
%   • totalValid(c,r): number of tested (stim → rec) pairs  
%   • totalSig(c,r):   number of those pairs with sigHigh == 1  
% Finally plots propMat with transparent NaNs.

% USER‐SPECIFIED PARAMETERS
subjDirs = {
    'subject 1'
    'subject 2'
    'subject 3'
    'subject 4'
    'subject 5'
    'subject 6'
    'subject 7'
};
binWidth =4;  % bin size in mm

% PASS 1: Gather all valid olf_coor across subjects for global bin edges
all_med_glob = [];
all_mid_glob = [];
for s = 1:numel(subjDirs)
    subjPath = subjDirs{s};
    if ~isfolder(subjPath), continue; end
    try
        load(fullfile(subjPath,'condensed_ecog_info_only.mat'),'condensed_ecog_info_only');
    catch, continue; end
    med = [condensed_ecog_info_only.lat_coor];
    mid = [condensed_ecog_info_only.mid_coor];
    olf = [condensed_ecog_info_only.olf_coor];
    valid = med(~isnan(med) & olf>-1);
    all_med_glob = [all_med_glob, valid];
end
if isempty(all_med_glob)
    error('No valid electrodes found across all subjects.');
end

minEdge = min(all_med_glob) - binWidth;
maxEdge = max(all_med_glob) + binWidth;

%overlapping bins
overlapFrac = 0;                        % fraction of the bin width that overlaps
step       = binWidth * (1 - overlapFrac);  
% define bin centers so that each bin is binWidth wide but starts step mm after the last
binCenters = (minEdge + binWidth/2) : step : (maxEdge - binWidth/2);
Nbins      = numel(binCenters);
% Create edges in steps of binWidth:
binEdges = minEdge : binWidth : maxEdge;


% Preallocate global accumulators
totalSig   = zeros(Nbins,Nbins);
totalValid = zeros(Nbins,Nbins);

% PASS 2: Loop over subjects to accumulate counts
for s = 1:numel(subjDirs)
    mirroredness_null=[];
    mirroredness=[];
    subjPath = subjDirs{s};
    if ~isfolder(subjPath)
        warning('Folder not found: %s. Skipping.', subjPath);
        continue;
    end
    % Load required data
    try
        load(fullfile(subjPath,'condensed_ecog_info_only.mat'),'condensed_ecog_info_only');
        load(fullfile(subjPath,'metadata.mat'),'metadata');
        load(fullfile(subjPath,'stimset.mat'),'stimset'); %#ok<NASGU>
    catch ME
        warning('Failed loading for %s: %s. Skipping.', subjPath,ME.message);
        continue;
    end
    % Gather valid electrodes
    med = [condensed_ecog_info_only.lat_coor];
    mid = [condensed_ecog_info_only.mid_coor];
    olf = [condensed_ecog_info_only.olf_coor];
    valid_idx = find(~isnan(med) & med~=0 & olf>-1);
    if isempty(valid_idx)
        warning('No valid electrodes in %s.', subjPath);
        continue;
    end
    med_valid = med(valid_idx);
    mid_valid = mid(valid_idx);
    % Bin assignment and side
    [~,~,binIdx] = histcounts(med_valid,binEdges);
    side_all = nan(size(mid_valid));
    side_all(mid_valid<0) = -1;
    side_all(mid_valid>0) = +1;

    % Build lookup
    electrodesInBin = cell(Nbins,2);
    for i = 1:numel(med_valid)
        d    = med_valid(i);
        side = side_all(i);
        % find all bin-centers whose window [c - bw/2, c + bw/2) contains d
        hits = find( d >= (binCenters - binWidth/2) & ...
            d <  (binCenters + binWidth/2) );
        for b = hits
            if side<0
                electrodesInBin{b,1}(end+1) = valid_idx(i);
            elseif side>0
                electrodesInBin{b,2}(end+1) = valid_idx(i);
            end
        end
    end

    % Compute subj counts
    subjSig   = zeros(Nbins,Nbins);
    subjValid = zeros(Nbins,Nbins);
    for c = 1:Nbins
        stimL = electrodesInBin{c,1};
        stimR = electrodesInBin{c,2};
        for r = 1:Nbins
            recL = electrodesInBin{r,1};
            recR = electrodesInBin{r,2};
            cntV = 0; cntS = 0;
            % Left→Right
            if ~isempty(stimL) && ~isempty(recR)
                for sElec = stimL
                    tg = metadata(sElec).targets;
                    if isempty(tg), continue; end
                    for rElec = recR
                        if rElec<=numel(tg) && isfield(tg(rElec),'sigHigh')
                            cntV = cntV+1;
                            cntS = cntS + tg(rElec).sigHigh;
                            mirroredness_null = [mirroredness_null;abs(med(rElec) - med(sElec))];
                            if tg(rElec).sigHigh
                                mirroredness = [mirroredness;abs(med(rElec) - med(sElec))];
                            else
                            end
                        end
                    end
                end
            end
            % Right→Left
            if ~isempty(stimR) && ~isempty(recL)
                for sElec = stimR
                    tg = metadata(sElec).targets;
                    if isempty(tg), continue; end
                    for rElec = recL
                        if rElec<=numel(tg) && isfield(tg(rElec),'sigHigh')
                            cntV = cntV+1;
                            cntS = cntS + tg(rElec).sigHigh;
                            mirroredness_null = [mirroredness_null;abs(med(rElec) - med(sElec))];
                            if tg(rElec).sigHigh
                                mirroredness = [mirroredness;abs(med(rElec) - med(sElec))];
                            else
                            end
                        end
                    end
                end
            end
            subjValid(r,c) = cntV;
            subjSig(r,c)   = cntS;
        end
    end
mirror_index(s) =  mean(mirroredness_null) - mean(mirroredness);

end

figure(1)
subplot(2,1,1)
histogram(mirroredness,'binwidth',5)
xlim([-50 50])
subplot(2,1,2)
histogram(mirroredness_null,'binwidth',5)
xlim([-50 50])
%%
figure(2)
hold on

% Scatter individual points with jitter
for i = 1:length(mirror_index)
    if ~isnan(mirror_index(i))
        plot(5 + rand, mirror_index(i), 'k.', 'MarkerSize', 20)
    end
end

% Compute summary stats
validVals = mirror_index(~isnan(mirror_index));
mVal = mean(validVals);
sdVal = std(validVals);

% Overlay mean as a red circle
plot(7, mVal, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8)

% Error bar for ± SD
errorbar(5.5, mVal, sdVal, 'r', 'LineWidth', 2, 'CapSize', 12)

% Axis formatting
ylim([-10 10])
xlim([2 8])
box off
grid off
set(gca, 'TickDir', 'out')

ylabel('Mirror Index')
title('Mirroredness Index per Subject')
