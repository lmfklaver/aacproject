%% Script to call for PETH RUN
% Long epochs only > 3 s
%
%
%
% % % % % % % % % % % % %
% % PETH RUN onset
% % % % % % % % % % % % %
figure
subplot(2,1,1)
timwin = [-5 5];
% Calculate
timeBefore  = abs(timwin(1));
timeAfter   = timwin(2);

selRunEpochs = runEpochs_long;

trlCenteredRunStart = selRunEpochs(:,1)-timeBefore;
trlCenteredRunStop  = selRunEpochs(:,1)+timeAfter;
trlCenteredRun      = [trlCenteredRunStart trlCenteredRunStop];

binSize     = options.binSize; % in sec
timeEdges   = timwin(1):binSize:timwin(2);
secsTot     = timwin(2)-timwin(1);

spike_toRun = realignSpikes(spikes, trlCenteredRun);

spikeTrl_Run = cell(1,length(selRunEpochs));

for iRun = 1:length(selRunEpochs)
    spikeTrl_Run{iRun} = spike_toRun{iAAC}{iRun} - selRunEpochs(iRun,1);
end

countHistoRun  = histcounts(cell2mat(spikeTrl_Run'),timeEdges);
rateHistoRun   = countHistoRun/length(selRunEpochs)*1/binSize; %
% maybe for time something like: linspace(timwin(1),timwin(2),length(rateHisto))



% Plot
histogram('BinEdges',timeEdges, 'BinCounts',rateHistoRun)
hold on
box off
set(gca,'TickDir','out')
title('PETH to Run')

xlabel('time(s)')
ylabel('spikes/s')
xlim(timwin)

%%
% % % % % % % % % % % % %
% % Raster Run
% % % % % % % % % % % % %
    subplot(2,1,2)

    plotSpkOffset = 0;

    for iRun = 1:length(selRunEpochs)
        selRunTr = spikeTrl_Run{iRun};
        plot(selRunTr',repmat(plotSpkOffset,1,length(selRunTr)),'k.');
        hold on
        plotSpkOffset = plotSpkOffset+1;
    end

    box off
    set(gca,'ydir','reverse')
    ylimits = get(gca,'YLim');
    xlabel('time (s)')
    ylabel('trials')
    set(gca,'TickDir','out')
    ylim([ylimits(1) plotSpkOffset])
    xlim(timwin)
