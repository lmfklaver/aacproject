% plot some spike width derivatives to show procedure

iAAC = 24;

figure, 
plot(spikes.filtWaveform{iAAC})

hold on
plot(diff(spikes.filtWaveform{iAAC}))
plot(diff(diff(spikes.filtWaveform{iAAC})))

plot(diff(diff(diff(spikes.filtWaveform{iAAC}))))
plot(diff(diff(diff(diff(spikes.filtWaveform{iAAC})))))

plot(diff(diff(diff(diff(diff(spikes.filtWaveform{iAAC}))))))

legend({'orginal','1st derivative','2nd derivative','3rd derivative', '4th derivative',' 5th derivative'})
%%
iAAC = 8;

figure, 
plot(spikes.filtWaveform{iAAC})

hold on
plot(diff(spikes.filtWaveform{iAAC}))
plot(diff(diff(spikes.filtWaveform{iAAC})))

plot(diff(diff(diff(spikes.filtWaveform{iAAC}))))
plot(diff(diff(diff(diff(spikes.filtWaveform{iAAC})))))
plot(diff(diff(diff(diff(diff(spikes.filtWaveform{iAAC}))))))

legend({'orginal','1st derivative','2nd derivative','3rd derivative', '4th derivative','5th derivative'})
title('filtered waveform')
