function [RipSpikePhase]=getRipSpikePhase(cd)
%%Modified code from Kaiser    

    ripspikes.times=[];
    basename = bz_BasenameFromBasepath(cd);
        load([basename '.ripples.events.mat'])
        load([basename '.gd_eps.mat'])
        load([basename '.spikes.cellinfo.mat'])
        [status,interval]=InIntervals(ripples.peaks(:,1),gd_eps); %Detect ripples outside of stim
        ripstart=ripples.timestamps(:,1);
        ripend=ripples.timestamps(:,2);
        gdrips=[];
        gdrips(:,1) = ripstart(status)-.05;
        gdrips(:,2) = ripend(status)+.05;
        [Congdrips] = ConsolidateIntervals(gdrips);
        for j = 1:length(spikes.times)
        [status] = InIntervals(spikes.times{j},Congdrips);
        ripspikes.times{j}=spikes.times{j}(status);
        end
    RipLFP = bz_GetLFP(ripples.detectorinfo.detectionparms.channel);
    filLFP = bz_Filter(RipLFP, 'passband', [120 250]);
    hilLFP = hilbert(double(filLFP.data));
    sigphaseLFP = angle(hilLFP);
    RipSpikePhase=[];
    for i=1:size(spikes.times,2)
        [N,BIN] = histc(ripspikes.times{i},RipLFP.timestamps);
        spiketimes=RipLFP.timestamps(logical(N));
        RipPhaseInd = ismember(RipLFP.timestamps,spiketimes);
        RipSpikePhase{i}=sigphaseLFP(RipPhaseInd);
    end
save([basename '.RipSpikePhase.analysis.mat'])
end
