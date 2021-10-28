%% Example theta trace + rasters AACs underneath
% 
% if strcmpi(basename,'u19_200313_155505')
%     opts.rippleChan     = 57;
% end

for iSess = [1:5,15:17]
    cd(dirN{iSess})
basepath = cd; basename = bz_BasenameFromBasepath(cd);
sessionInfo = bz_getSessionInfo;
load([basename '.ripples.events.mat'])

rippleChan = ripples.detectorinfo.detectionparms.channel;
data = bz_LoadBinary([basename '.dat'],...
    'frequency',30000,'nChannels', sessionInfo.nChannels,'channels',...
    rippleChan+1,'downsample',30');
% lfp = bz_GetLFP(rippleChan);
% spikes = bz_LoadPhy;

%%
selXLim = [3010 3020];
figure,
set(gcf,'PaperOrientation','landscape');
set(gcf,'Position',[50 50 1200 400]);

% s1 = subplot(3,1,1);
% plot(lfp.timestamps,lfp.data)
% xlim(selXLim)

% get and plot theta
% lfp_theta = BandpassFilter(double(lfp.data), 1250 ,[1 100]); % nb changed from [6 8]
lfp_theta= double(data);
time = (1:length(data))/1000;
% plot(lfp.timestamps,lfp_theta*0.195)
p1 = plot(time, lfp_theta*0.195);
xlim([selXLim])


for iThEpoch = 1:100:time(end)
    widthTE = 1;
    thEpoch = [iThEpoch-widthTE iThEpoch+widthTE];
    xlim(thEpoch);
    
    unitStr = ['theta_example_' num2str(iThEpoch) '.pdf'];
    print(gcf,unitStr, '-dpdf','-bestfit')
    append_pdfs(['E:\Dropbox\PD_Hpc\Progress\' basename '_Raw_20201216.pdf'],[unitStr])
    delete([unitStr])
end


end
