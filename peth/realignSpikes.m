function [spikes_realigned] = realignSpikes(spikes, trl)
% inputs: 
% spikes struct from buzcode
% trl = [Nx2] start stop of trials 
% example: [spike_toPulse] = realignSpikes(spikes, trlCenteredEpoch);


spikes_realigned = cell(1,length(spikes.times));

for iUnit  = 1:length(spikes.times)
    [status, interval]= InIntervals(spikes.times{iUnit},trl);
    numSpkAligned = histoc(interval(interval>0),1:size(trl,1));
    spikes_realigned{iUnit} = mat2cell(spikes.times{iUnit}(status),numSpkAligned);
end

end
