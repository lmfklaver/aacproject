%%Find spikes in intervals
%% goodeps/not good eps - softcode pulsetime-loop through sessions

noRunEpochs = noRunTimes;
noPulseEpochs = noPulseTimes;
totalLengthStill = sum(noRunEpochs(:,2)-noRunEpochs(:,1));

% runEpochs = [2200 2500];
%
totalLengthRun = sum(runEpochs(:,2)-runEpochs(:,1));

%run
[status, ~] = cellfun(@(a) InIntervals(a, runEpochs),spikes.times,'uni', false);
for iUnit = 1:length(spikes.times)
 runSpikes{iUnit} = spikes.times{iUnit}(status{iUnit});
end 
[status, ~] = cellfun(@(a) InIntervals(a, pulseEpochs),runSpikes,'uni', false);
for iUnit = 1:length(spikes.times)
 runSpikesPulse{iUnit} = spikes.times{iUnit}(status{iUnit});
end
[status, ~] = cellfun(@(a) InIntervals(a, noPulseEpochs),runSpikes,'uni', false);
for iUnit = 1:length(spikes.times)
 runSpikesNoPulse{iUnit} = spikes.times{iUnit}(status{iUnit});
end

%noRun
[status, ~] = cellfun(@(a) InIntervals(a, noRunEpochs),spikes.times,'uni', false);
for iUnit = 1:length(spikes.times)
 noRunSpikes{iUnit} = spikes.times{iUnit}(status{iUnit});
end
[status, ~] = cellfun(@(a) InIntervals(a, pulseEpochs),noRunSpikes,'uni', false);
for iUnit = 1:length(spikes.times)
 noRunSpikesPulse{iUnit} = spikes.times{iUnit}(status{iUnit});
end
[status, ~] = cellfun(@(a) InIntervals(a, noPulseEpochs),noRunSpikes,'uni', false);
for iUnit = 1:length(spikes.times)
 noRunSpikesNoPulse{iUnit} = spikes.times{iUnit}(status{iUnit});
end


% %total pulses time with spikes in it? 
% % runFR all
% for iUnit = 1:length(spikes.times)
% fr_allrun(iUnit) = length(runSpikes{iUnit})/totalLengthRun;
% end
% 
% mFRallrun = mean(fr_allrun(~INTIndx));
% stdFRallrun = std(fr_allrun(~INTIndx));
%%

[pyrs, ints, aacs] = splitCellTypes(basepath)

%runFR pulse
totalPulseTimeRec = length(pulseEpochs)*0.3;%300ms time here
for iUnit = 1:length(spikes.times)
fr_allrunpulse(iUnit) = length(runSpikesPulse{iUnit})/totalPulseTimeRec; %firing rate within pulse total number of timestamps/total number of timestamps, entire pulse is bin
end
mFRallrunpulse = mean(fr_allrunpulse(pyrs) );
stdFRallrunpulse = std(fr_allrunpulse(pyrs)./sqrt(length(spikes.times)));

%runFR nopulse
totalTimeNoPulseRec =totalLengthRun - totalPulseTimeRec;
for iUnit = 1:length(spikes.times)
fr_allrunnopulse(iUnit) = length(runSpikesNoPulse{iUnit})/totalTimeNoPulseRec;
end
mFRallrunnopulse = mean(fr_allrunnopulse(pyrs));
stdFRallrunnopulse = std(fr_allrunnopulse(pyrs)./sqrt(length(spikes.times)));

%norunFR pulse
totalTimeNoPulseRec =totalLengthRun - totalPulseTimeRec;
for iUnit = 1:length(spikes.times)
fr_allnorunpulse(iUnit) = length(noRunSpikesPulse{iUnit})/totalPulseTimeRec;
end
mFRallnorunpulse = mean(fr_allnorunpulse(pyrs));
stdFRallnorunpulse = std(fr_allnorunpulse(pyrs)./sqrt(length(spikes.times)));



%norunFR nopulse
totalTimeNoPulseRec =totalLengthRun - totalPulseTimeRec;
for iUnit = 1:length(spikes.times)
fr_allnorunnopulse(iUnit) = length(noRunSpikesNoPulse{iUnit})/totalTimeNoPulseRec;
end
mFRallnorunnopulse = mean(fr_allnorunnopulse(pyrs));
stdFRallnorunnopulse = std(fr_allnorunnopulse(pyrs)./sqrt(length(spikes.times)));


%% Total Time RECORDINGS is now not corrected for how many pulses fell during running and during rest!! 
%% NEEDS to BE CORRECTED!!! 



%%
figure
bar([mFRallrunpulse,mFRallnorunpulse,mFRallrunnopulse,mFRallnorunnopulse])
hold on
e1=errorbar([mFRallrunpulse,mFRallnorunpulse,mFRallrunnopulse,mFRallnorunnopulse],...
    [[stdFRallrunpulse,stdFRallnorunpulse,stdFRallrunnopulse,stdFRallnorunnopulse]]);
set(e1,'LineStyle','none')
set(gca,'XTickLabel',{'run pulse','norun pulse','run nopulse','norun nopulse'})
% 
% [status, ~] = cellfun(@(a) InIntervals(a, stillEpochs),spikes.times,'uni', false);
% for iUnit = 1:length(spikes.times)
%  stillSpikes{iUnit} = spikes.times{iUnit}(status{iUnit});
% end
%%
signrank(fr_allrunpulse(~INTIndx),fr_allnorunpulse(~INTIndx))
signrank(fr_allrunnopulse(~INTIndx),fr_allnorunnopulse(~INTIndx))
signrank(fr_allrunpulse(~INTIndx),fr_allrunnopulse(~INTIndx))
signrank(fr_allnorunpulse(~INTIndx),fr_allnorunnopulse(~INTIndx))

%%
gain_norun =(fr_allnorunpulse(~INTIndx)-fr_allnorunnopulse(~INTIndx))./fr_allnorunnopulse(~INTIndx);
gain_run = (fr_allrunpulse(~INTIndx)-fr_allrunnopulse(~INTIndx))./fr_allrunnopulse(~INTIndx);

edges = -1:0.1:2;
norun_n = histcounts(gain_norun,edges);
run_n = histcounts(gain_run,edges);

sumP = sum(norun_n);
sumQ = sum(run_n);

figure,

normalizedq = run_n./sumQ*100;
bar(normalizedq)
hold on
normalizedp = norun_n./sumP*100;
bar(normalizedp) % not in pulse
hold on
box off
xlabel('Gain')
ylabel('% Occurance')
title('Distribution of Gain PYR''s RUN v NO RUN')
legend({'RUN','NO RUN'})
xt = get(gca,'XTick');
xl = edges(xt(2:end));
set(gca,'XTick',xt(2:end),'XTickLabel', string(xl))



print(gcf,'distributionRUNnoRUN.pdf','-dpdf')

%%
% plot ([gain]
figure,
plot([gain_run;gain_norun],'ko-')
xlim([0.5 2.5])
ylim([-1.5 3])
xt = get(gca,'XTick');
set(gca,'XTick',[1 2],'XTickLabel', {'RUN','NO RUN'})
box off
ylabel('Gain')
