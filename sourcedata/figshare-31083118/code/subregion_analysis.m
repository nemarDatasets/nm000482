subjDirs = {
    'subject 1'
    'subject 2'
    'subject 3'
    'subject 4'
    'subject 5'
    'subject 6'
    'subject 7'
};


Ncomp = 4;
colors = lines(Ncomp);

overall_proportion_sig = zeros(Ncomp,Ncomp);

for comp = 1:Ncomp
    allCoords = [];
    allSig = [];
    allLats = [];
    subject = [];
    proportion_sig_subj = zeros(length(subjDirs),Ncomp);

    for s = 1:numel(subjDirs)
        subjPath = subjDirs{s};
        if ~isfolder(subjPath), continue; end

        try
            load(fullfile(subjPath,'condensed_ecog_info_only.mat'),'condensed_ecog_info_only');
            load(fullfile(subjPath,'metadata.mat'),'metadata');
        catch
            continue;
        end

        olf = [condensed_ecog_info_only.olf_coor];
        valid_idx = find(olf > -1);

        for sElec = valid_idx
            stimComp = getCompartment(sElec,condensed_ecog_info_only);
            stimSide = condensed_ecog_info_only(sElec).mid_coor > 0;
            if stimComp ~= comp, continue; end

            tg = metadata(sElec).targets;
            if isempty(tg), continue; end

            for rElec = 1:numel(tg)
                if rElec > numel(condensed_ecog_info_only), continue; end

                recComp = getCompartment(rElec,condensed_ecog_info_only);
                recSide = condensed_ecog_info_only(rElec).mid_coor > 0;

                % Only cross-hemispheric
                if recSide == stimSide, continue; end

                isSig = isfield(tg(rElec), 'sigHigh') && tg(rElec).sigHigh;
                allSig(end+1)   = isSig;
                allCoords(end+1)= recComp;
                allLats(end+1)  = tg(rElec).onset_latency;
                subject(end+1) = s;
            end
        end
    end

    % Preallocate results
    proportion_sig = nan(1,4);
    avg_lats       = nan(1,4);
    lats_error     = nan(1,4);

    for compartment = 1:4
        proportion_sig(compartment) = ...
            sum(allSig==1 & allCoords==compartment) / ...
            sum(allCoords==compartment);

        avg_lats(compartment) = ...
            mean(allLats(allSig==1 & allCoords==compartment),'omitnan');

        lats_error(compartment) = ...
            2 * std(allLats(allSig==1 & allCoords==compartment), 'omitnan') / ...
            sqrt(sum(allSig==1 & allCoords==compartment));

        for t = 1:length(subjDirs)
            proportion_sig_subj(t,compartment) = sum(allSig==1 & allCoords==compartment  & subject == t) / ...
            sum(allCoords==compartment & subject == t);
        end

        if compartment == comp
            figure(compartment)
lat_data = allLats(allSig==1 & allCoords==compartment);
[mean(lat_data) std(lat_data)]

% Create a tiled layout for stacked plots with shared X
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% --- Histogram on top ---
nexttile
histogram(lat_data, 'BinWidth', 3)
xlim([8 50])
ylim([0 10])
box off
grid off
set(gca, 'TickDir', 'out')
ylabel('Count')
title(sprintf('Compartment %d', compartment))

% --- Boxplot on bottom (horizontal) ---
nexttile
boxplot(lat_data, 'orientation', 'horizontal', 'Notch', 'on', 'Widths', 0.5)
xlim([8 50])
set(gca, 'TickDir', 'out', 'YTick', [])
xlabel('Latency (ms)')
lat_data_all{comp}=lat_data;

        else
        end
    end



    overall_proportion_sig(comp,:) = proportion_sig;

end
%%
figure(2)
imagesc(overall_proportion_sig)
clim([0 0.6])

nC = 256;
whiteToRed = [ linspace(1,0,nC)'*0.5 + 0.5 , ...  % red channel
               linspace(1,0,nC)', ...             % green channel
               linspace(1,0,nC)' ];               % blue channel
colormap(whiteToRed);

box off
grid off
set(gca, 'TickDir', 'out')
colorbar  

%%


function c = getCompartment(elec, info)
c = 5;  % default to unclassified

if info(elec).med_coor<0 && info(elec).olf_coor>-1
    c = 1;  % Medial OFC
elseif info(elec).med_coor >0 &&  info(elec).lat_coor <0 && info(elec).trans_coor>0
    c = 2;  % Anterior OFC
elseif info(elec).med_coor >0 &&  info(elec).lat_coor <0 && info(elec).trans_coor<0
    c = 3;  % Posterior OFC
elseif info(elec).lat_coor>0
    c = 4;  % Lateral OFC
elseif info(elec).olf_coor<0
        c = 0;  % Ventromedial OFC
    end
end



