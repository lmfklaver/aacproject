% conert pulseEpochs to pulses events mat and optoSim manipulation mat for
% compatibility with CellExplorer
% make sure to check default opsin etc. 


pulses.timestamps  = pulseEpochs;
pulses.peaks            = (pulseEpochs(:,2)-pulseEpochs(:,1))/2;
pulses.amplitude        = [];
pulses.amplitudeUnits   = [];
pulses.eventID          = ones(length(pulseEpochs),1);
pulses.eventIDlabels    = cell(1,length(pulses.timestamps));
pulses.eventIDlabels(:) = {'OptoStim - Arch'};
% pulses.eventIDlabels: cell array with labels for classifying various event types defined in stimID (cell array, Px1).
% pulses.eventIDbinary: boolean specifying if eventID should be read as binary values (default: false).
pulses.center           = pulses.peaks;%;
pulses.duration         = pulseEpochs(:,2)-pulseEpochs(:,1);
pulses.detectorinfo     = 'getPulseEpochs'; 

save(strcat(basename, '.pulses.events.mat'), 'pulses')

optoStim.timestamps = pulses.timestamps;
optoStim.stimID     = pulses.eventID;%(==2)
optoStim.center     = pulses.center;
optoStim.duration   = pulses.duration;

save([basename '.optoStim.manipulation.mat'],'optoStim')