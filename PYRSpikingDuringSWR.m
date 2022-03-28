%% plot of number of Pyramidal Cells in SWR stim vs control PYR/AACS

load([basename '.ripples.events.mat']);
load([basename '.optoStim.manipulation.mat']);
load([basename '.spikes.cellinfo.mat']);
load([basename '_celltypes.mat']);
pulseEpochs = optoStim.timestamps;

[peakInPulse, pulseWithRip] = findRipplesInPulse(ripples, pulseEpochs);


pyrLog = cellfun(@(a) contains(a, 'pyr'),allcelltypes);
pyrInd = find(pyrLog);
%%
for iPyr = pyrInd
    [status,interval] = InIntervals(spikes.times{iPyr},ripples.timestamps(pulseWithRip,:));
    v = 1:length(ripples.timestamps);
    RipOutsidePulse = ~ismember(v,pulseWithRip);
    
        [statusNO,intervalNO] = InIntervals(spikes.times{iPyr},ripples.timestamps(RipOutsidePulse,:));

    numSpkperRip_ON(iPyr) = sum(status);
    numSpkperRip_OFF(iPyr) = sum(statusNO);
    
    % per ripple
    
        for iInterval = unique(interval(interval~=0))';
            numSpkperRip_ONper(iPyr,iInterval) = sum(ismember(iInterval,interval));
        end
        for iIntervalNO = unique(intervalNO(intervalNO~=0))';
            numSpkperRip_OFFper(iPyr,iIntervalNO) = sum(ismember(iIntervalNO,intervalNO));
        end
    % gemiddelde spikes per ripple, opslaan numSpks per rip?
    
end

compareONOFF = [numSpkperRip_ON',numSpkperRip_OFF'];

scatter(numSpkperRip_ON',numSpkperRip_OFF')
xlim([0 40000])
%%

figure,plot(repmat([1,2],length(numSpkperRip_ON)),[numSpkperRip_ON/length(numSpkperRip_ON)',numSpkperRip_OFF/length(numSpkperRip_OFF)],'ko-')
xlim([0 3])
title('NumSpikes during SWRs Excitation AACs 1 Session')
set(gca,'XTick',[1,2],'XTickLabel',{'ON','OFF'})
ylabel(['Average # of Spikes In SWR' num2str(length(spikes.times)) ' PYR cells'])
%%
figure
 histogram(compareONOFF(:,1)),hold on, histogram(compareONOFF(:,2));