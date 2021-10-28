lfp = bz_GetLFP('all');
xml = LoadXml(fullfile([basename,'.xml']));

timwin = [-0.5 0.5];
timeBefore = abs(timwin(1));
timeAfter = timwin(2);

trlCenteredPulseStart = pulseEpochs(:,1)-timeBefore;
trlCenteredPulseStop = pulseEpochs(:,1)+timeAfter;

trlCenteredPulse = [trlCenteredPulseStart trlCenteredPulseStop];

data_uVo = double(lfp.data)*0.0195; % to microvolts; %samples x channels
% lpfilter < 250Hz
data_uV = lowpass(data_uVo,250,1250);
%%
for iAnatGrps = 1:length(xml.AnatGrps)
    %%
    figure
    set(gcf,'Position',[680    85   560   893])
    
    plotOffset =0;
    for iChan = 1:length(xml.AnatGrps(iAnatGrps).Channels)
        
        selChannel = xml.AnatGrps(iAnatGrps).Channels(iChan)+1; %or 0idx
         selData = data_uV(:,selChannel);
        
        [status, intervals] = InIntervals(lfp.timestamps, trlCenteredPulse);
        
       for iPulse = 1:length(unique(intervals))-1
           selDataPulse = selData(intervals==iPulse);
        lfp_pulseTrials(iPulse,:) = selDataPulse(1:1250); %soms 1251 door downsample mismatch in ts
       end

       if selChannel == 5
                  plotOffset = plotOffset - 20;

           continue
           
       else
       plot(mean(lfp_pulseTrials)+ plotOffset,'k'),
       hold on
       plotOffset = plotOffset - 20;
       end
    end
    % stim period rectangle
    box off
    
    r1 = rectangle('Position',[625 20 375 5],'FaceColor',[0 0 0]);
    t1 = text(625,40, 'Stim On');
    r2 = rectangle('Position',[1250 -305 10 10],'FaceColor',[0 0 0]);
    t2 = text(1265,-300,'10 uV');
    
    ax1 = gca;
ax1.YAxis.Visible = 'off';

title(['pulse-triggered LFPs of shank ' num2str(iAnatGrps)])
xlabel(['time in samples' ])
% % time = linspace(-0.5,0.5,1251);
% % timSize = size(lfp_pulseTrials,2)+1;
% % xt      = 1:250:1250;
% % xl      =time(xt)
% % xl = 
% % strxl   = string(xl);
% % 
% % set(gca,'XTick', xt,'XTickLabel',strxl)
     %%
    print(gcf,['pulseTriggeredLFP_shank' num2str(iAnatGrps)],'-dpdf')
    %time 
    
end
       
        