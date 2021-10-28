function [baserate, runrate] = getRatesTrialsBaseRun(spikes, runEpochs, timwin,options)

binSize     = options.binSize; % in sec
%%
runwidth = runEpochs(1,2)-runEpochs(1,1);

[spkTimesPerTrial] = getSpkTimTrials(spikes, runEpochs, timwin,options); %0 is start run

for iUnit = 1:length(spikes.UID)
    baseInd{iUnit} = cellfun(@(a) a<0, spkTimesPerTrial{iUnit},'UniformOutput', false);
    stimInd{iUnit} = cellfun(@(a) a>0 & a <runwidth, spkTimesPerTrial{iUnit},'UniformOutput', false);
end

for iUnit = 1:length(spikes.UID)
    for iTr = 1:length(runEpochs)
        
        baseIndTr = baseInd{iUnit}{iTr};
        stimIndTr = stimInd{iUnit}{iTr};
        
        selUnitTr = spkTimesPerTrial{iUnit}{iTr};
        selUnitBaseTr = selUnitTr(baseIndTr);
        selUnitStimTr = selUnitTr(stimIndTr);
        
        basetrials{iUnit}{iTr}=selUnitBaseTr;
        stimtrials{iUnit}{iTr}=selUnitStimTr;
        
        baserate_temp{iUnit}{iTr} = length(basetrials{iUnit}{iTr})*1/abs(timwin(1));
        %assuming that timwin1 is negative to 0 , baseline calculated over entire timwin negative
        stimrate_temp{iUnit}{iTr} = length(stimtrials{iUnit}{iTr})*1/runwidth;%
    end
end

for iUnit = 1:length(spikes.UID)
    baserate{iUnit} = cell2mat(baserate_temp{iUnit})';
    runrate{iUnit} = cell2mat(stimrate_temp{iUnit})';
end
   
end
    