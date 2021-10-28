% ChR Sessions

% folderAACchr = 'E:\Dropbox\Unc5B\Spks_ait_rip_for_lianne';
% folderAACchr = 'E:\Data\AxoAxo\Spks_ait_rip_for_lianne' % this one misses
% files

folderAACchr = 'E:\Dropbox\PD_Hpc\Unc5B\Spks_ait_rip_for_lianne';
sessions = {'mouse1_180414','mouse1_180415','mouse1_180501a','mouse1_180501b',...
    'mouse1_180502a','mouse1_180502b','mouse3_180627'...
    'mouse3_180628','mouse3_180629','mouse4_181114b',...
    'mouse5_181112B','mouse5_181116','mouse6_190330', ...
    'mouse6_190331'};
 cd(folderAACchr)
% opts.binSize=0.01
%%
rate_concat = [];
resc_rate_concat = [];
for iSess = 1:length(sessions)
    selecSession = sessions{iSess};
    
    
    
    [pulseEpochs, ripples] = getPulseAndRipsFromEvt(selecSession);
    pulseDur = pulseEpochs(:,2)-pulseEpochs(:,1);
    pulseEpochs(pulseDur<0.09,:) = []; % to get those weird ones out 
    
% % %     cd('D:\Data\Analogin')
% % %     analogin_file   = [sessions{iSess}, '_analogin.dat'];
% % %     basename = selecSession;
% % %     [analogin.pulse, analogin.pos, analogin.reward, analogin.ts] = getAnaloginVals(basename,params,board_adc_channels,opts);
% % %     
% % %     %% Get Pulse Epochs
% % %     [pulseEpochs] = getPulseTimes(analogin);
% % %     
    
    % RipplePeak in Pulse?
    [peakInPulse, pulseWithRip, aRP] = inh_getPeakInPulse(ripples.peaks, pulseEpochs);
    
    % Load Spikes
    load([folderAACchr filesep selecSession '.spikes.cellinfo.mat']);
    %getCCG
    [ccg,t]=CCG(spikes.times,[],'Fs',params.sampFreq, 'binSize',opts.ccgBinSize,'duration', opts.ccgDur, 'norm', 'rate');
    
    timwin = [-1 1];
    % all pulses
    [rate, ~, ~] = inh_rastersToPulse_new(spikes, pulseEpochs,timwin,ccg,t,opts);
    rowmin = min(rate,[],2);
    rowmax = max(rate, [],2);
    
    resc_rate = rescale(rate,'InputMin',rowmin, 'InputMax', rowmax);
    
    rate_concat = [rate_concat;rate];
    resc_rate_concat = [resc_rate_concat;resc_rate];
    
% %     AACIdx_concatenated = []; % make criteria here
end
%%
Z = zscore(rate_concat,[],2);
Zs = sortrows(Z, 1010,'ascend');
figure,imagesc(Zs)
xlim([900 1100])
caxis([0 5])
box off

InclInd_FR = mean(rate_concat,2)<3;
InclInd_Inh = mean(rate_concat(:,1000:1100),2)<mean(rate_concat(:,900:1100),2); % rate stim smaller than rate pulse
IncInd = InclInd_FR & InclInd_Inh;
figure,imagesc(Z(IncInd,:))
% inh_heatmap_rates(rate_concat, windowLength,INTIdx, AACIdx)
