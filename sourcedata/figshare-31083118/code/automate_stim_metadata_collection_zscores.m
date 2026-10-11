% List of subject folders (relative to current directory)
% this goes through each trace and stores relevant metadata
% (whether there are significant responses in response to stim given a particular stimulated electrode and recorded electrode)
% run this before other codes
subjDirs = {
    'subject 1'
    'subject 2'
    'subject 3'
    'subject 4'
    'subject 5'
    'subject 6'
    'subject 7'
};

zThresh = 5;

% Define window and alpha once
window = [200, 500];  % [ms before, ms after]

for iSubj = 1:numel(subjDirs)
    subjPath = subjDirs{iSubj};
    
    % Verify folder exists
    if ~isfolder(subjPath)
        warning('Folder not found: %s. Skipping.', subjPath);
        continue;
    end
    
    % Attempt to load condensed_ecog and stimset from this subject folder
    try
        dataCE = load(fullfile(subjPath, 'condensed_ecog_OFC_only.mat'), 'condensed_ecog');
        dataSS = load(fullfile(subjPath, 'stimset.mat'), 'stimset');
    catch ME
        warning('Failed loading necessary files for %s: %s. Skipping.', subjPath, ME.message);
        continue;
    end
    
    condensed_ecog = dataCE.condensed_ecog;
    stimset        = dataSS.stimset;
    
    % Initialize metadata struct array for this subject
    nCh = numel(condensed_ecog);
    metadata = repmat(struct('targets', []), nCh, 1);
    
    % Loop over each stimulating electrode
    for eStim = 1:length(stimset.electrode)
        
        if isempty(stimset.electrode(eStim).intensity)
            continue;
        end
        
        % --- find highest intensity index ---
        nInt = numel(stimset.electrode(eStim).intensity);
        ints = nan(1, nInt);
        for k = 1:nInt
            ints(k) = stimset.electrode(eStim).intensity(k).amplitude;
        end
        [~, idxMax] = max(ints);
        timesHigh = stimset.electrode(eStim).intensity(idxMax).times;
        nHigh     = numel(timesHigh);
        if nHigh == 0
            continue;
        end
        
        % --- Build PSTH for the highest intensity ---
        nTime = window(1) + window(2) + 1;
        psthHigh = repmat(struct('data', []), nCh, 1);
        for eRec = 1:nCh
            if any(isnan(condensed_ecog(eRec).bipolar))
                continue;
            end
            psthMat = zeros(nHigh, nTime);
            for t = 1:nHigh
                tt = timesHigh(t);
                baselineIdx = (tt - window(1)):(tt -1);
                respIdx     = (tt - window(1)):(tt + window(2));
                bl = mean(condensed_ecog(eRec).bipolar(baselineIdx));
                segment = condensed_ecog(eRec).bipolar(respIdx) - bl;
                psthMat(t, :) = segment;
            end
            psthHigh(eRec).data = psthMat;
        end
        
        % --- Determine zero index and searchStart (>10 ms) ---
        zeroIdx    = find((-window(1):window(2)) == 0);  % = window(1)+1
        searchStart = zeroIdx +10;
        searchEnd = zeroIdx + 50;
        
        % -----------------------------
        % SIGNIFICANCE & LATENCY LOGIC
        % -----------------------------
        sigHigh     = false(nCh,1);
        onsetLatency = nan(nCh,1);
        peakLatency = nan(nCh,1);
        peakAmp = nan(nCh,1);
        
        for eRec = 1:nCh
            if isempty(psthHigh(eRec).data)
                continue;
            end
            
            dataMat   = psthHigh(eRec).data;      % [nHigh × nTime]
            meanTrace = mean(dataMat, 1);         % 1×nTime
            
            % Compute mean trace and z-score
            base = meanTrace(1:zeroIdx);
            mu   = mean(base);
            sd   = std(base);
            zt   = (meanTrace - mu) / sd;

            % Detect significance
            idxRel = find(zt(searchStart:searchEnd) > zThresh, 1, 'first');
            
            if ~isempty(idxRel)
                sigHigh(eRec) = true;
                onsetIdx = searchStart + idxRel - 1;
                onsetLatency(eRec) = onsetIdx - zeroIdx;
                [peakAmp(eRec), pkRel] = max(zt(searchStart:searchEnd));
                peakIdx = searchStart + pkRel - 1;
                peakLatency(eRec) = peakIdx - zeroIdx;
            end
        end


        % -----------------------------
        % POPULATE METADATA FOR THIS eStim
        % -----------------------------
        metadata(eStim).targets = repmat(struct( ...
            'peakAmp',   [], ...
            'onset_latency',   [], ...
            'peak_latency',   [], ...
            'sigHigh',   [] ), nCh, 1);
        
        for eRec = 1:nCh
            metadata(eStim).targets(eRec).peakAmp   = peakAmp(eRec, :);
            metadata(eStim).targets(eRec).onset_latency   = onsetLatency(eRec);
            metadata(eStim).targets(eRec).peak_latency   = peakLatency(eRec);
            metadata(eStim).targets(eRec).sigHigh   = sigHigh(eRec);
        end


    end
    
    % Save the computed metadata.mat inside this subject's folder
    save(fullfile(subjPath, 'metadata.mat'), 'metadata', '-v7.3');
end