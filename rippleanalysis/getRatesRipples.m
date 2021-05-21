function [ratepertrial] = getRatesRipples(spikes, ripples, options)


trlCenteredRipStart = ripples.timestamps(:,1);
trlCenteredRipStop = ripples.timestamps(:,2);

trlCenteredRip = [trlCenteredRipStart trlCenteredRipStop];

%%
binSize     = options.binSize; % in sec

%%
spike_toPulse = realignSpikes(spikes, trlCenteredRip);

%%
for iUnit = 1:length(spikes.UID)
    plotSpkOffset = 0;
    
    for iPulse = 1:length(ripEpochs)
        spikeTrl_Pulse{iUnit}{iPulse} = spike_toPulse{iUnit}{iPulse} - ripEpochs(iPulse,1);
    end
end
ratepertrial = spikeTrl_Pulse;

end
    