% Arch Sessions

% folderAACchr = 'E:\Dropbox\Unc5B\Spks_ait_rip_for_lianne';
% folderAACchr = 'E:\Data\AxoAxo\Spks_ait_rip_for_lianne' % this one misses
% files

folderAACunc5b = 'D:\Data\Axoaxonic_Data_Lianne';
sessions = {'u19_200310_135409','u19_200313_120452','u19_200313_155505',...
    'u21_200305_153604','m175_200821_151859'};

cd(folderAACchr)
% opts.binSize=0.001
%%
rate_concatArch = [];
resc_rate_concatArch = [];
for iSess = 1:length(sessions)
    selecSession = sessions{iSess};
    folderAACsel = [folderAACunc5b filesep selecSession];
    cd(folderAACsel)
    
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
%     [peakInPulse, pulseWithRip, aRP] = inh_getPeakInPulse(ripples.peaks, pulseEpochs);
    
    % Load Spikes
    load([folderAACsel filesep selecSession '.spikes.cellinfo.mat']);
    %getCCG
    [ccg,t]=CCG(spikes.times,[],'Fs',params.sampFreq, 'binSize',opts.ccgBinSize,'duration', opts.ccgDur, 'norm', 'rate');
    
    timwin = [-1 1];
    % all pulses
    [rate, ~, ~] = inh_rastersToPulse_new(spikes, pulseEpochs,timwin,ccg,t,opts);
    rowmin = min(rate,[],2);
    rowmax = max(rate, [],2);
    
    resc_rate = rescale(rate,'InputMin',rowmin, 'InputMax', rowmax);
    
    rate_concatArch = [rate_concatArch;rate];
    resc_rate_concatArch = [resc_rate_concatArch;resc_rate];
    
% %     AACIdx_concatenated = []; % make criteria here
end
%%
Z = zscore(rate_concatArch,[],2);
Zs = sortrows(Z, 1010,'ascend');
figure,imagesc(Zs)
xlim([900 1100])
caxis([0 5])
box off

InclInd_FR = mean(rate_concatArch,2)<3;
InclInd_Inh = mean(rate_concatArch(:,1000:1100),2)<mean(rate_concatArch(:,900:1100),2); % rate stim smaller than rate pulse
IncInd = InclInd_FR & InclInd_Inh;
figure,imagesc(Z(IncInd,:))
% inh_heatmap_rates(rate_concat, windowLength,INTIdx, AACIdx)
