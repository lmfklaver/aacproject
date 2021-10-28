function [ratepertrial] = getRatesTrials(spikes, pulseEpochs, timwin,options)

timeBefore = abs(timwin(1));
timeAfter = timwin(2);

trlCenteredPulseStart = pulseEpochs(:,1)-timeBefore;
trlCenteredPulseStop = pulseEpochs(:,1)+timeAfter;

trlCenteredPulse = [trlCenteredPulseStart trlCenteredPulseStop];

%%
binSize     = options.binSize; % in sec
timeEdges   = timwin(1):binSize:timwin(2);
secsTot     = timwin(2)-timwin(1);
%     secsPerBin = secsTot/numBins; % for when hist in stead of histcounts

%%
spike_toPulse = realignSpikes(spikes, trlCenteredPulse);

%%
for iUnit = 1:length(spikes.UID)
    plotSpkOffset = 0;
    
    for iPulse = 1:length(pulseEpochs)
        spikeTrl_Pulse{iUnit}{iPulse} = spike_toPulse{iUnit}{iPulse} - pulseEpochs(iPulse,1);
    end
end
ratepertrial = spikeTrl_Pulse;




% 
%         baseSpikes              = spikeTrl_Pulse{iPulse}(spikeTrl_Pulse{iPulse}<0);
%         pulseTimeIndices        = find(spikeTrl_Pulse{iPulse}>0 & spikeTrl_Pulse{iPulse}<.300);
%         pulseSpikes             = spikeTrl_Pulse{iPulse}(pulseTimeIndices);
%         mean_bl_rate{iUnit}{iPulse}     = length(baseSpikes)/abs(timwin(1));
%         puldur_rate{iPulse}     = length(pulseSpikes)/.3;
%         
% %         ratechange_pulse{iPulse} = puldur_rate./mean_bl_rate; % niet in procenten
%     end
%     
%     [P(iUnit),H(iUnit)] = signrank(cell2mat(mean_bl_rate),cell2mat(puldur_rate));
%     
%     countHisto = histcounts(cell2mat(spikeTrl_Pulse'),timeEdges);
%     rateHisto   = countHisto/length(pulseEpochs)*1/binSize; %
%     timeHisto   = 1:length(rateHisto); % fix this still
%     % maybe for time something like: linspace(timwin(1),timwin(2),length(rateHisto))
%     
%     rate(iUnit,:)   = rateHisto;
%     count(iUnit,:)  = countHisto;
%     time(iUnit,:)   = timeHisto;
%     ratesBLStim  = [];    
end
    