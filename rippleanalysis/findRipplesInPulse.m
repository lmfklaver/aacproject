function [peakInPulse, pulseWithRip] = findRipplesInPulse(ripples, pulseEpochs)

% turn into function: findRipsInPulse
% 2/28/2022 ripple entirely in pulse

pulseWithRip = [];
peakInPulse = [];
for iPeak = 1:length(ripples.peaks)
    ripInPulse = find(ripples.timestamps(iPeak,1)>=pulseEpochs(:,1) & ripples.timestamps(iPeak,2)<=pulseEpochs(:,2));
    if ~isempty(ripInPulse)
        pulseWithRip = [pulseWithRip ripInPulse];
        peakInPulse = [peakInPulse iPeak];
    end
end
end

