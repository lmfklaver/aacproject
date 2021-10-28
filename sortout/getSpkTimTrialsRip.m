function [spkTimesPerTrial] = getSpkTimTrialsRip(spikes, ripples)
% This function
%
%   USAGE 
%
%   Dependencies: 
%
%   INPUTS
%   Name-value paired inputs:
%  
%   OUTPUTS
%   
%   EXAMPLES
%
%
%   NOTES
%  
%
%   TO-DO
%  
%
%   HISTORY
% 

trlRip = ripples.timestamps;

spike_toRip = realignSpikes(spikes, trlRip);

%%
for iUnit = 1:length(spikes.UID)
    plotSpkOffset = 0;
    
    for iRip = 1:length(trlRip)
        spikeTrl_Rip{iUnit}{iRip} = spike_toRip{iUnit}{iRip} - trlRip(iRip,1);
    end
end

spkTimesPerTrial = spikeTrl_Rip;


end
    