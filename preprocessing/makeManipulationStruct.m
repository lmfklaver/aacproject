%% makes pulses from manipulationStruct

% load([basename '.pulses.events.mat'])
% pulseEpochs = optoStim.timestamps;

pulses.timestamps  = pulseEpochs;
pulses.peaks = (pulseEpochs(:,2)-pulseEpochs(:,1))/2;
pulses.amplitude = [];
pulses.amplitudeUnits = [];
pulses.eventID = ones(length(pulseEpochs),1);
pulses.eventIDlabels = cell(1,length(pulses.timestamps));
pulses.eventIDlabels(:) = {'OptoStim - ChR'};
% pulses.eventIDlabels: cell array with labels for classifying various event types defined in stimID (cell array, Px1).
% pulses.eventIDbinary: boolean specifying if eventID should be read as binary values (default: false).
pulses.center = pulses.peaks;%;
pulses.duration = pulseEpochs(:,2)-pulseEpochs(:,1);
pulses.detectorinfo = 'getPulseEpochs'; 

save(strcat(basename, '.pulses.events.mat'), 'pulses')

% % for excluding manipulations in cell explorer: make manipulation.mat

optoStim.timestamps = pulses.timestamps;
optoStim.stimID = pulses.eventID;%(==2)
optoStim.center = pulses.center;
optoStim.duration = pulses.duration;


save([basename '.optoStim.manipulation.mat'],'optoStim')

% find(cell_metrics.optoStim_modulationIndex>1.5)

session = sessionTemplate(basepath,'showGUI',true);
cell_metrics = ProcessCellMetrics('session', session);