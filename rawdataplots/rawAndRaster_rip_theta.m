%% Example theta trace + rasters AACs underneath
% 
% if strcmpi(basename,'u19_200313_155505')
%     opts.rippleChan     = 57;
% end
basepath = cd; basename = bz_BasenameFromBasepath(cd);

load([basename '.ripples.events.mat'])

rippleChan = ripples.detectorinfo.detectionparms.channel;
lfp = bz_GetLFP(rippleChan);
spikes = bz_LoadPhy;

%%
selXLim = [3010 3020];
figure,
% s1 = subplot(3,1,1);
% plot(lfp.timestamps,lfp.data)
% xlim(selXLim)

% get and plot theta
s2 = subplot(2,1,1);
lfp_theta = BandpassFilter(double(lfp.data), 1250 ,[1 100]); % nb changed from [6 8]
plot(lfp.timestamps,lfp_theta*0.195)
xlim([selXLim])
xlabel('time (s)')
ylabel('amplitude (uV)')

% now plote AAC spike times in it
s3 =subplot(2,1,2);
% startUnder = -800%min(lfp_theta);
% plotShift = 0;
[~, ~, AACs] = splitCellTypes(basepath);

for iSel = 1:length(AACs)
    iAAC = AACs(iSel);
    x = spikes.times{iAAC};
    totSpkTims = length(spikes.times{iAAC});
    y = [repmat(iSel-1,totSpkTims,1),repmat(iSel,totSpkTims,1)];
    hold on
    
    tx = [x.';x.';nan(1,length(x))];
    ty = [y(:,1).';y(:,2).';nan(1,length(x))];
    plot(tx(:),ty(:))
    xlim([selXLim])
    xlabel('time (s)')
    ylabel('AAC')
end

for iThEpoch = 1:100:lfp.timestamps(end)
    widthTE = 1;
    thEpoch = [iThEpoch-widthTE iThEpoch+widthTE];
    s1.XLim=thEpoch;
    s2.XLim=thEpoch;
    s3.XLim=thEpoch;
    
            unitStr = ['theta_example_' num2str(iThEpoch) '.pdf'];
    print(gcf,unitStr, '-dpdf','-bestfit')
        append_pdfs(['E:\Dropbox\PD_Hpc\Progress\1_100_Raw_Session_5_20201215.pdf'],[unitStr])
        delete([unitStr])
end


%%
for iRip = 1:length(ripples.peaks)
    
    selRipple = iRip;
    timeAroundripple = .2;
    ripXLim = [ripples.peaks(selRipple)-timeAroundripple ripples.peaks(selRipple)+timeAroundripple];
    
    figure,
    %       subplot(3,1,1)
    %       plot(lfp.timestamps,lfp.data)
    %       xlim(ripXLim)
    
    % get and plot theta
    subplot(2,1,1)
    % get and plot ripple
    lfp_ripple = BandpassFilter(double(lfp.data), 1250 ,[100 250]);
    plot(lfp.timestamps,lfp_ripple*0.195)
    
    xlim(ripXLim)
    xlabel('time s')
    ylabel('amplitude uV')
    subplot(2,1,2)
    
    for iSel = 1:length(AACs)
        iAAC = AACs(iSel);
        x = spikes.times{iAAC};
        totSpkTims = length(spikes.times{iAAC});
        y = [repmat(iSel-1,totSpkTims,1),repmat(iSel,totSpkTims,1)];
        hold on
        
        tx = [x.';x.';nan(1,length(x))];
        ty = [y(:,1).';y(:,2).';nan(1,length(x))];
        plot(tx(:),ty(:))
        xlim(ripXLim)
        xlabel('time s')
        ylabel('AAC')
    end
    
    print(gcf,['ripple_example_' num2str(selRipple) '.pdf'], '-dpdf','-bestfit')
    close
end