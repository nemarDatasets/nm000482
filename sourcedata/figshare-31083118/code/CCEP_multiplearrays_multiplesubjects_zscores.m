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
binWidth = 5;  % bin size in mm

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
overlapFrac = 0.5;                        % fraction of the bin width that overlaps
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
                            if tg(rElec).sigHigh
                            [sElec rElec]
                            else
                            end
                        end
                    end
                end
            end
            %Right→Left
            if ~isempty(stimR) && ~isempty(recL)
                for sElec = stimR
                    tg = metadata(sElec).targets;
                    if isempty(tg), continue; end
                    for rElec = recL
                        if rElec<=numel(tg) && isfield(tg(rElec),'sigHigh')
                            cntV = cntV+1;
                            cntS = cntS + tg(rElec).sigHigh;
                            if tg(rElec).sigHigh
                            [sElec rElec]
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
    % Accumulate
    totalValid = totalValid + subjValid;
    totalSig   = totalSig   + subjSig;
end

% FINAL: Compute and plot with opacity proportional to data density
propMat = totalSig ./ totalValid;
propMat(totalValid==0) = NaN;


% Number of colors in the map
nC = 256;

% build white→red colormap:
%   R = ones(nC,1)
%   G = linspace(1,0,nC)'
%   B = linspace(1,0,nC)'



set(gca,'YDir','normal');
% custom white→red
nC = 256;
whiteToRed = [ linspace(1,0,nC)'*0.5 + 0.5 , linspace(1,0,nC)',linspace(1,0,nC)'];
colormap(whiteToRed);

set(gca,'Color',[1 1 1]);
alphaMask = totalValid / max(totalValid(:));  % 0 = transparent, 1 = opaque
hImg = imagesc(binCenters, binCenters, propMat);
%hImg.AlphaData = alphaMask;
colorbar;
caxis([0 1]);

xticks([-30:10:60]);
yticks([-30:10:60]);
xlim([-35 30])
ylim([-35 30])
xtickangle(45);

xlabel('Stimulus‐bin center (mm)');
ylabel('Recording‐bin center (mm)');
%%
% Mask valid tiles
validMask = totalValid > 0 & ~isnan(propMat);

% Observed significant counts
obsSig = round(propMat(validMask) .* totalValid(validMask));

% Expected counts under null = proportional to totalValid
totalSigCount = sum(obsSig);
expSig = totalSigCount * ( totalValid(validMask) / sum(totalValid(validMask)) );

% Chi-square test
chi2stat = sum((obsSig - expSig).^2 ./ expSig);
df = numel(obsSig) - 1;
pVal = 1 - chi2cdf(chi2stat, df);

fprintf('Chi-square = %.2f, df = %d, p = %.4g\n', chi2stat, df, pVal);

%%
% Find valid tiles
validMask = totalValid > 0 & ~isnan(propMat);

% Row and column indices of valid tiles
[rowIdx, colIdx] = find(validMask);
distances = abs(rowIdx - colIdx);  % distance from diagonal

% ----- Actual significant connections -----
sigDistances = [];
for i = 1:numel(rowIdx)
    r = rowIdx(i); c = colIdx(i);
    nSig = round(propMat(r,c) * totalValid(r,c)); % number sig
    sigDistances = [sigDistances; repmat(distances(i), nSig, 1)];
end

% ----- Null distribution -----
nullDistances = [];
for i = 1:numel(rowIdx)
    nAll = totalValid(rowIdx(i), colIdx(i));
    nullDistances = [nullDistances; repmat(distances(i), nAll, 1)];
end

% ----- Compare -----
[pVal,~,stats] = ranksum(sigDistances, nullDistances); % Wilcoxon rank-sum
fprintf('Median dist (sig) = %.2f, Median dist (null) = %.2f, p = %.4g\n', ...
    median(sigDistances), median(nullDistances), pVal);

% Optional: plot histograms
figure;
histogram(nullDistances, 'FaceColor', [0.7 0.7 0.7], 'DisplayName', 'Null', 'BinWidth', 1); hold on
histogram(sigDistances, 'FaceColor', [0.9 0.2 0.2], 'DisplayName', 'Significant', 'BinWidth', 1);
xlabel('Distance from diagonal (bins)');
ylabel('Count');
legend;
box off;
set(gca,'TickDir','out');




%% LATENCIES

clear binEdges binIdx_all Nbins electrodesInBin propMat binCenters binLabels
clear minEdge maxEdge side_all
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
binWidth = 5;  % bin size in mm

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
    olf = [condensed_ecog_info_only.olf_coor];
    valid = med(~isnan(med) & med~=0 & olf>-1);
    all_med_glob = [all_med_glob, valid];
end
if isempty(all_med_glob)
    error('No valid electrodes found across all subjects.');
end

minEdge = min(all_med_glob) - binWidth;
maxEdge = max(all_med_glob) + binWidth;

%overlapping bins
overlapFrac = 0.5;                        % fraction of the bin width that overlaps
step       = binWidth * (1 - overlapFrac);  
% define bin centers so that each bin is binWidth wide but starts step mm after the last
binCenters = (minEdge + binWidth/2) : step : (maxEdge - binWidth/2);
Nbins      = numel(binCenters);
% Create edges in steps of binWidth:
binEdges = minEdge : binWidth : maxEdge;


% Preallocate global accumulators
totalSig   = zeros(Nbins,Nbins);
totalValid = zeros(Nbins,Nbins);
totalLat = zeros(Nbins,Nbins);

% PASS 2: Loop over subjects to accumulate counts
for s = 1:numel(subjDirs)
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
    subjAvglat = zeros(Nbins,Nbins);
    for c = 1:Nbins
        stimL = electrodesInBin{c,1};
        stimR = electrodesInBin{c,2};
        for r = 1:Nbins
            recL = electrodesInBin{r,1};
            recR = electrodesInBin{r,2};
            cntV = 0; cntS = 0;latS = 0;
            % Left→Right
            if ~isempty(stimL) && ~isempty(recR)
                for sElec = stimL
                    tg = metadata(sElec).targets;
                    if isempty(tg), continue; end
                    for rElec = recR
                        if rElec<=numel(tg) && isfield(tg(rElec),'sigHigh')
                            cntV = cntV+1;
                            cntS = cntS + tg(rElec).sigHigh;
                            if tg(rElec).sigHigh
                                latS = latS + tg(rElec).onset_latency;else;end
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
                            if tg(rElec).sigHigh
                                tg(rElec).onset_latency
                                latS = latS + tg(rElec).onset_latency; else;end
                        end
                    end
                end
            end
            subjValid(r,c) = cntV;
            subjSig(r,c) = cntS;
            subjAvglat(r,c)   = latS;
        end
    end
    % Accumulate
    totalValid = totalValid + subjValid;
    totalSig   = totalSig   + subjSig;
    totalLat = totalLat + subjAvglat;
end

% FINAL: Compute and plot with opacity proportional to data density
latMat = totalLat ./ totalSig;
latMat(totalSig==0) = NaN;

% Create alpha mask (opacity) scaled by totalValid
alphaMask = totalSig/max(totalSig(:))>0;  % 0 = transparent, 1 = opaque

figure;
% Plot propMat and capture image handle
% draw the image so that its columns map to binCenters on X,
% and its rows to binCenters on Y:
hImg = imagesc(binCenters, binCenters, latMat);
set(gca,'YDir','normal');      % put row 1 at the bottom
colormap(flipud(jet));
hImg.AlphaData = alphaMask;
set(gca,'Color',[1 1 1]);
colorbar;
caxis([10 50]);
axis normal;                   

% now ticks are automatically at your binCenters,
% and you can fine‐tune them if you like:
xticks(binCenters);
yticks(binCenters);
xtickangle(45);

xlabel('Stimulus‐bin center (mm)');
ylabel('Recording‐bin center (mm)');


%% Latency distribution (interhemispheric pairs only)

subjDirs = {
    'subject 1'
    'subject 2'
    'subject 3'
    'subject 4'
    'subject 5'
    'subject 6'
    'subject 7'
};



allOnsets = [];

for s = 1:numel(subjDirs)
    p = subjDirs{s};
    if ~isfolder(p), continue; end

    try
        M = load(fullfile(p,'metadata.mat'),'metadata');
        C = load(fullfile(p,'condensed_ecog_info_only.mat'),'condensed_ecog_info_only');
    catch
        continue;
    end

    meta = M.metadata;
    info = C.condensed_ecog_info_only;

    for eStim = 1:numel(meta)

        % Need stim coordinate to define hemisphere
        if eStim > numel(info), continue; end
        stimMid = info(eStim).mid_coor;
        if isnan(stimMid) || stimMid==0, continue; end
        stimSide = stimMid > 0;   % 0=left (mid<0), 1=right (mid>0)

        tg = meta(eStim).targets;
        if isempty(tg), continue; end

        for eRec = 1:numel(tg)

            % Need rec coordinate to define hemisphere
            if eRec > numel(info), continue; end
            recMid = info(eRec).mid_coor;
            if isnan(recMid) || recMid==0, continue; end
            recSide = recMid > 0;

            % Only cross-hemispheric
            if recSide == stimSide, continue; end

            if isfield(tg(eRec),'sigHigh') && tg(eRec).sigHigh && isfield(tg(eRec),'onset_latency')
                allOnsets(end+1) = tg(eRec).onset_latency; %#ok<AGROW>
            end
        end
    end
end

figure;
histogram(allOnsets, 'BinWidth', 2);   % actually 2-ms bins
xlabel('Onset latency (ms)');
ylabel('Count');
title('Distribution of CCEP Onset Latencies (Contralateral Only)');



