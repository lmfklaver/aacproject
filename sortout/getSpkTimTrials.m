function [spkTimesPerTrial] = getSpkTimTrials(spikes, pulseEpochs, timwin,options)

timeBefore = abs(timwin(1));
timeAfter = timwin(2);

trlCenteredPulseStart = pulseEpochs(:,1)-timeBefore;
trlCenteredPulseStop = pulseEpochs(:,1)+timeAfter;

trlCenteredPulse = [trlCenteredPulseStart trlCenteredPulseStop];

%%
binSize     = options.binSize; % in sec

%%
spike_toPulse = realignSpikes(spikes, trlCenteredPulse);

%%
for iUnit = 1:length(spikes.UID)
    plotSpkOffset = 0;
    
    for iPulse = 1:length(pulseEpochs)
        spikeTrl_Pulse{iUnit}{iPulse} = spike_toPulse{iUnit}{iPulse} - pulseEpochs(iPulse,1);
    end
end
spkTimesPerTrial = spikeTrl_Pulse;


end
    