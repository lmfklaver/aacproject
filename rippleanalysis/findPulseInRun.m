function [pulseinrun] = findPulseInRun(runepochs, pulseEpochs)

% Findpulseinrun 2022

runwithpulse = [];
pulseinrun = [];
for iPeak = 1:length(pulseEpochs)
    ripInPulse = find(pulseEpochs(iPeak,1)>=runepochs.run.epochs(:,1) & pulseEpochs(iPeak,2)<=runepochs.run.epochs(:,2));
    if ~isempty(ripInPulse)
        runwithpulse = [runwithpulse ripInPulse];
        pulseinrun = [pulseinrun iPeak];
    end
end
end

