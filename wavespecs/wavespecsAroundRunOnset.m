%% Wavespecs around run onset
%
%
%   Use Runepochs and Wavespec to calculate average wavespecs around run
%   onset


load([basename '_analogin.mat'])
load([basename '.wavespec.analysis.mat'])

[vel] = getVelocity(analogin, 'downsampleFactor',300,'circDisk',...
2*pi*26,'doFigure',true);

% check normalizations, why does on one wheel the animal appears to run 10
% times faster than on the other?

minvel_cm_s = 10;

[run] = getRunEpochs(basepath,vel,'minRunSpeed',minvel_cm_s);

%%
xmln    = [basepath filesep basename '.xml'];
fname   = [basepath filesep basename '.lfp'];
xml     = LoadXml(xmln);

load([basename '.ripples.events.mat'])
ch      = ripples.detectorinfo.detectionchannel;
lfp     = bz_GetLFP(ch);

%% cut snippets of wavespec that correspond to runEpochs out 
selRunEpochs = run.epochs; % in seconds
selRunIdx = run.index; % matching indices

timMS   = 1; % s
timSamp = timMS * lfp.samplingRate;
ops.tw_ws =  timSamp;
% ops.bl_ws =  ops.tw_ws*2;

ws_time =[];
ws_data = [];
countRuns = 0;
for iRun = 1:length(selRunEpochs)
    selRunStartIdx = selRunIdx(iRun,1);
    %nb this is to exclude first or last epoch that might fall outside the
    %window of interest)
    if abs(selRunStartIdx-ops.tw_ws) == selRunStartIdx-ops.tw_ws
        countRuns = countRuns + 1;
        %       if abs(lfp.timestamps(selRunStart)-ops.bl_ws) == lfp.timestamps(selRunStart)-ops.bl_ws
        ws_time(:,:,countRuns) = wavespec.timestamps(selRunStartIdx-ops.tw_ws:selRunStartIdx+ops.tw_ws,:);
        ws_data(:,:,countRuns) = wavespec.data(selRunStartIdx-ops.tw_ws:selRunStartIdx+ops.tw_ws,:);
    end
    
end


wsd_m = mean(ws_data,3);

figure,
imagesc(abs(wsd_m'))
set(gca,'YDir','normal')

if length(wsd_m)/2 ~=round(length(wsd_m)/2)
    tsamp_mid = (length(wsd_m)+1)/2;
else
    tsamp_mid = (length(wsd_m))/2;
end

tsamp_start = 1;
tsamp_stop = length(wsd_m);

xt = [tsamp_start tsamp_mid tsamp_stop];
strxl = {['-' num2str(timMS)], 0, num2str(timMS)}
set(gca,'XTick', xt,'XTickLabel',strxl)

ylabel('frequency')
xlabel('time(s)')
