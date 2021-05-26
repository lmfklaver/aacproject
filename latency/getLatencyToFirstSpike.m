function LatencyFirstSpike = getLatencyToFirstSpike(spikeTrl_Pulse, pulseEpochs)
% This function is meant to calculate the latency to the First spike after
% optostim
%
%   USAGE
%   LatencyFirstSpike = getLatencyToFirstSpike(spikeTrl_Pulse, pulseEpochs)
%
%   Dependencies:
%   
%
%   INPUTS
%   spkTrl_Pulse    - spiketimes centered to pulse is 0, struct with
%                      {iUnit}{iPulse}
%   pulseEpochs     - Nx2 matrix of start and stop times 
%
%
%   OUTPUTS
%     LatencyFirstSpike - Nx1 For each trial in PulseEpochs the latency that it
%     took for the first spike to occur
%
%   EXAMPLES
%   LatencyFirstSpike = getLatencyToFirstSpike(peth.trials{iUnit},pulseEpochs)
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