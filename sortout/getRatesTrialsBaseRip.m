function [riprate] = getRatesTrialsBaseRip(spikes, ripples)
% This function
%
%   USAGE 
%
%   Dependencies: 
%   Buzcode, Englishlab\utilities
%
%   INPUTS
%   spikes
%   ripples
%
%   Name-value paired inputs:
%   
%   OUTPUTS
%   riprate
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
%   2021/1  Lianne is documenting this thing 

%%

ripLength = ripples.timestamps(:,2) - ripples.timestamps(:,1); 
spkTimesPerTrial = getSpkTimTrialsRip(spikes, ripples);

% spikes not in ripple and not in pulse


% for iUnit = 1:length(spikes.UID)
%     ripInd{iUnit} = cellfun(@(a) a>0 & a <.150, spkTimesPerTrial{iUnit},'UniformOutput', false);
% end

for iUnit = 1:length(spikes.UID)
    for iRip = 1:length(ripLength)
        
%         stimIndTr = ripInd{iUnit}{iRip};
        
        selUnitTr = spkTimesPerTrial{iUnit}{iRip};
%         selUnitStimTr = selUnitTr(stimIndTr);
        
        riptrials{iUnit}{iRip}=selUnitTr;
        riprate_temp{iUnit}{iRip} = length(riptrials{iUnit}{iRip})*1/ripLength(iRip);%
    end
end

for iUnit = 1:length(spikes.UID)
%     baserate{iUnit} = cell2mat(baserate_temp{iUnit})';
    riprate{iUnit} = cell2mat(riprate_temp{iUnit})';
end
   
end
    