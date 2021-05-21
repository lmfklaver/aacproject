function LatencyFirstSpike = getLatencyToFirstSpike(spikeTrl_Pulse, pulseEpochs)
% This function is meant to calculate the latency to the First spike after
% optostim
%
%   USAGE
%   
%
%   Dependencies:
%   Buzcode
%
%   INPUTS
%   Name-value paired inputs:
%
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
%   HISTORY
%   2021/1  Lianne made this into a function - was a script before
%
%

%% Get Latency


for iPulse = 1:length(pulseEpochs)
    positiveSpk = spikeTrl_Pulse{iPulse}(spikeTrl_Pulse{iPulse}>0);
    
    if ~isempty(positiveSpk)
        LatencyFirstSpike(iPulse) = positiveSpk(1);
    else
        LatencyFirstSpike(iPulse) = NaN; % for an empty trialß
    end
end
end