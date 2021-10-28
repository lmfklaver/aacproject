% Get ripple characteristics inside and outside pulse
% e.g. amp, freq, etc.


sessions =   {'D:\Data\Axoaxo\u19_200310_135409';...%16
            'D:\Data\Axoaxo\u19_200313_155505'};...%17

ripplespNPin = []; ripplespNPout = [];
datapFin = []; datapFout=[];
datapAin = []; datapAout = [];
dataDurin = []; dataDurout = [];
%%
for iSess = 2%1:length(sessions)

%%
    basepath = sessions{iSess};
    cd(basepath)
    basename = bz_BasenameFromBasepath(basepath);
    load([basename '.ripples.events.mat']);
    load([basename '.optoStim.manipulation.mat']);
    
    %
    intervals   = optoStim.timestamps;
    
    [statusPeak,~]          = InIntervals(ripples.peaks,intervals);
    [statusStart,intStart]  = InIntervals(ripples.timestamps(:,1),intervals);
    [statusStop,~]          = InIntervals(ripples.timestamps(:,2),intervals);

    delayRipToPulse = ripples.timestamps(statusStart,1)-intervals(intStart(statusStart),1);
    figure,
    histogram(delayRipToPulse)
    % pick one 
%     status = statusPeak;
    status = statusStart & statusStop;
 
    ripchan     = ripples.detectorinfo.detectionchannel;
    lfp         = bz_GetLFP(ripchan);
    filtered    = bz_Filter(double(lfp.data),'filter','butter','passband',[100 250],'order', 3);
    timestamps  = lfp.timestamps;
    
    [maps,data,stats] = bz_RippleStats(filtered,timestamps,ripples);
    
    ripplespNPin    = [ripplespNPin; ripples.peakNormedPower(status)];
    ripplespNPout   = [ripplespNPout; ripples.peakNormedPower(~status)];
    datapFin        = [datapFin; data.peakFrequency(status)];
    datapFout       = [datapFout; data.peakFrequency(~status)];
    datapAin        = [datapAin; data.peakAmplitude(status)];
    datapAout       = [datapAout;data.peakAmplitude(~status)];
    dataDurin       = [dataDurin;data.duration(status)];
    dataDurout      = [dataDurout;data.duration(~status)];
end
%%
load([basename '.spikes.cellinfo.mat'])
    load([basename '.optoStim.manipulation.mat']);

for iUnit = 1:length(spikes.times)
        [statusSpk,~]          = InIntervals(spikes.times{iUnit},optoStim.timestamps);
        spksInPulse(iUnit) = length(spikes.times{iUnit}(statusSpk));
        spksOutPulse(iUnit) = length(spikes.times{iUnit}(~statusSpk));
end


figure,
subplot(4,1,1)
histogram(ripplespNPin)
hold on
histogram(ripplespNPout)
title('ripples.peakNormedPower')
legend({'Rip in Pulse','Rip out Pulse'})

% how this is calculated, from bz_FindRipples
%     signal = bz_Filter(double(p.Results.lfp),'filter','butter','passband',passband,'order', 3);
% squaredSignal = signal.^2
% [normalizedSquaredSignal,sd] = unity(Filter0(window,sum(squaredSignal,2)),sd,keep);
% max of normalizedSquaredSignal is what is the peakNormPower 
% U = (A - meanA)/stdA;
% max Value is added 


subplot(4,1,2)
histogram(datapFin)
hold on 
histogram(datapFout)
title('data.peakFrequency')
legend({'Rip in Pulse','Rip out Pulse'})

subplot(4,1,3)
histogram(datapAin)
hold on 
histogram(datapAout)
title('data.peakAmplitude')
legend({'Rip in Pulse','Rip out Pulse'})
%data.peakAmplitude = maps.amplitude(:,centerBin);
% maps.amplitude = SyncMap(a,i,'durations',durations,'nbins',nBins,'smooth',0);

subplot(4,1,4)
histogram(dataDurin)
hold on 
histogram(dataDurout)
title('data.duration')
legend({'Rip in Pulse','Rip out Pulse'})
%%

figure
subplot(4,2,1)
imagesc(maps.phase(status,:))
title('instantaneous phase (one ripple per row)')
subplot(4,2,2)
imagesc(maps.phase(~status,:))
title('instantaneous phase (one ripple per row)')
subplot(4,2,3)
imagesc(maps.amplitude(status,:))
title('envelope amplitude (one ripple per row)')
subplot(4,2,4)
imagesc(maps.amplitude(~status,:))
title('envelope amplitude (one ripple per row)')
subplot(4,2,5)
imagesc(maps.ripples(status,:))
title('instantaneous amplitude (one ripple per row)')
subplot(4,2,6)
imagesc(maps.ripples(~status,:))
title('instantaneous amplitude (one ripple per row)')
subplot(4,2,7)
imagesc(maps.frequency(status,:))
title('instantaneous frequency (one ripple per row)')
subplot(4,2,8)
imagesc(maps.frequency(~status,:))
title('instantaneous frequency (one ripple per row)')
%%
figure,
subplot(2,1,1)
plot(stats.acg.t,stats.acg.data)
subplot(2,1,2)
