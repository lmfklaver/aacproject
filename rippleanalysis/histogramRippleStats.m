% Get ripple characteristics inside and outside pulse
% e.g. amp, freq, etc.

basepath = cd;
basename = bz_BasenameFromBasepath(basepath);
load([basename '.ripples.events.mat']);
load([basename '.optoStim.manipulation.mat']);

%
values = ripples.peaks;
intervals = optoStim.timestamps;

[status,interval] = InIntervals(values,intervals);

%%

%% 
ripchan = ripples.detectorinfo.detectionchannel;

lfp = bz_GetLFP(ripchan);
filtered = BandpassFilter(double(lfp.data), 1250,[100 250]);
timestamps = lfp.timestamps;

[maps,data,stats] = bz_RippleStats(filtered,timestamps,ripples);

%%
figure,
subplot(4,1,1)
histogram(ripples.peakNormedPower(status))
hold on
histogram(ripples.peakNormedPower(~status))
title('ripples.peakNormedPower')
legend({'Rip in Pulse','Rip out Pulse'})

subplot(4,1,2)
histogram(data.peakFrequency(status))
hold on 
histogram(data.peakFrequency(~status))
title('data.peakFrequency')
legend({'Rip in Pulse','Rip out Pulse'})

subplot(4,1,3)
histogram(data.peakAmplitude(status))
hold on 
histogram(data.peakAmplitude(~status))
title('data.peakAmplitude')
legend({'Rip in Pulse','Rip out Pulse'})

subplot(4,1,4)
histogram(data.duration(status))
hold on 
histogram(data.duration(~status))
title('data.duration')
legend({'Rip in Pulse','Rip out Pulse'})
%%

figure
subplot(4,1,1)
imagesc(maps.phase)
title('instantaneous phase (one ripple per row)')
subplot(4,1,2)
imagesc(maps.amplitude)
title('envelope amplitude (one ripple per row)')
subplot(4,1,3)
imagesc(maps.ripples)
title('instantaneous amplitude (one ripple per row)')
subplot(4,1,4)
imagesc(maps.frequency)
title('instantaneous frequency (one ripple per row)')
%%
figure,
plot(stats.acg.t,stats.acg.data)
