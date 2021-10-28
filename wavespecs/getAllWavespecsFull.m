function [out] = getAllWavespecsFull(dirN)

% Calc all wavespecs for sessions for run check
for iSess = [1:5, 8, 9, 15:length(dirN)]
    
    cd(dirN{iSess})
    basepath = cd; basename = bz_BasenameFromBasepath(cd);
    
    load([basename '.ripples.events.mat'])
    rippleChan      = ripples.detectorinfo.detectionchannel;
    
    lfp = bz_GetLFP(rippleChan);

    freqRange = [1 200];
    numFreqs = freqRange(end)-freqRange(1)*2;
    
    wavespec         = bz_WaveSpec(lfp,'frange',freqRange,'nfreqs',numFreqs,'space','lin');
%     wavespec.data    = abs(wavespec.data);
    ws_reshaped     = reshape(abs(wavespec.data),[length(wavespec.timestamps),wavespec.nfreqs,countRuns]);
    
    save([basename '.wavespec.analysis.mat'],'wavespec','ws_reshaped');
    out = 1;
end
