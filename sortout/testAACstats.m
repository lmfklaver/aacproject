%% Test AACs out


options.binSize = 0.001;
% 30000 samplespersec 

timwin = [-.5 5];
pulsewidth = pulseEpochs(1,2)-pulseEpochs(1,1);

[spkTimesPerTrial] = getRatesTrials(spikes, pulseEpochs, timwin,options);

for iUnit = 1:length(spikes.UID)
    baseInd{iUnit} = cellfun(@(a) a<0, spkTimesPerTrial{iUnit},'UniformOutput', false);
    stimInd{iUnit} = cellfun(@(a) a>0 & a <pulsewidth, spkTimesPerTrial{iUnit},'UniformOutput', false);
end

% spikes per base or stim

for iUnit = 1:length(spikes.UID)
for iTr = 1:length(pulseEpochs)
    
    baseIndTr = baseInd{iUnit}{iTr};
    stimIndTr = stimInd{iUnit}{iTr};
    
    selUnitTr = spkTimesPerTrial{iUnit}{iTr};
    selUnitBaseTr = selUnitTr(baseIndTr);
    selUnitStimTr = selUnitTr(stimIndTr);
    
    basetrials{iUnit}{iTr}=selUnitBaseTr;
    stimtrials{iUnit}{iTr}=selUnitStimTr;
    
    baserate{iUnit}{iTr} = length(basetrials{iUnit}{iTr})/abs(timwin(1));
    stimrate{iUnit}{iTr} = length(stimtrials{iUnit}{iTr})/abs(pulsewidth);
end
end

% a2 = cellfun(@(x) length(x),a);

for iUnit = 1:length(spikes.UID)

baseCount(iUnit,:) = cell2mat(cellfun(@(a) length(a), basetrials{iUnit},'uni',false'));
stimCount(iUnit,:) = cell2mat(cellfun(@(a) length(a), stimtrials{iUnit},'uni',false'));

baseHz(iUnit,:) = cell2mat(baserate{iUnit});
stimHz(iUnit,:) = cell2mat(stimrate{iUnit});

end

%simple tttest
for iUnit = 1:length(spikes.UID),[h(iUnit) p(iUnit)] = ttest(baseCount(iUnit,:),stimCount(iUnit,:));,end
ttestp= find(p<0.05)

%signrank
for iUnit = 1:length(spikes.UID),[p(iUnit) h(iUnit)] = signrank(baseHz(iUnit,:),stimHz(iUnit,:));,end
signrankp =find(p<0.05)
for iUnit = 1:length(spikes.UID),[p(iUnit) h(iUnit)] = signrank(baseHz(iUnit,:),stimHz(iUnit,:)),'tail','right');,end
signranktailp =find(p<0.05)

% signtest
for iUnit = 1:length(spikes.UID),[p(iUnit) h(iUnit)] = signtest(baseHz(iUnit,:),stimHz(iUnit,:));,end
n = find(p<0.05)


% 
% % stats
% for iUnit = 1:length(spikes.UID)
%     
% allTimesBase = cell2mat(baserate{iUnit}(:));
% allTimesStim = cell2mat(stimrate{iUnit}(:));
% 
% [h(iUnit), p(iUnit)] = ttest2(allTimesStim, allTimesBase,'Tail','right');
% 
% % signrank
% 
% end

% other things: Zeta, Salt

