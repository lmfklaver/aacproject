% for SALT

load([basename '.spikes.cellinfo.mat'])
% load([basename '.pulses.events.mat'])
load([basename '.optoStim.manipulation.mat'])
pulseEpochs = optoStim.timestamps;

options.binSize = 1/30000;%0.001;
% 30000 samplespersec 
%%
timwin = [-.5 5];
pulsewidth = pulseEpochs(1,2)-pulseEpochs(1,1);


[spkTimesPerTrial] = getRatesTrials(spikes, pulseEpochs, timwin,options);

for iUnit = 1:length(spikes.UID)
    baseInd{iUnit} = cellfun(@(a) a>-0.48 & a<0, spkTimesPerTrial{iUnit},'UniformOutput', false);
    stimInd{iUnit} = cellfun(@(a) a>0 & a <0.01 , spkTimesPerTrial{iUnit},'UniformOutput', false);
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
    end
end



% binarize

for iUnit = 1:length(spikes.UID)
    for iTr = 1:length(pulseEpochs)    
        rasteredgesbase = timwin(1):options.binSize:0;
        rasteredgesstim = 0:options.binSize:pulsewidth;
        binaryBase{iUnit}(iTr,:) = histcounts(basetrials{iUnit}{iTr},rasteredgesbase);
        binaryStim{iUnit}(iTr,:) = histcounts(stimtrials{iUnit}{iTr},rasteredgesstim);
    end
end

for iUnit = 1:length(spikes.UID)
spt_baseline = binaryBase{iUnit};
spt_test = binaryStim{iUnit};
dt = options.binSize;


[p(iUnit) I(iUnit)] = salt(spt_baseline,spt_test,dt,0.001);
end
find(p<0.01)




