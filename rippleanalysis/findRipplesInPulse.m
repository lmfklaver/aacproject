function [peakInPulse, pulseWithRip] = findRipplesInPulse(ripples, pulseEpochs)

% turn into function: findRipsInPulse


pulseWithRip = [];
peakInPulse = [];

for iPeak = 1:length(ripples.peaks)
    ripInPulse = find(ripples.peaks(iPeak)>pulseEpochs(:,1) & ripples.peaks(iPeak)<pulseEpochs(:,2)) ;
    if ~isempty(ripInPulse)
        pulseWithRip = [pulseWithRip ripInPulse];
        peakInPulse = [peakInPulse iPeak];
    end
end
end
