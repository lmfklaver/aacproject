rmpath(genpath('E:\Dropbox\Code\Plotting'))

dirN  = {...
    'D:/Data/AxoAxo/mouse1/mouse1_180412_2';... %2
    'D:/Data/AxoAxo/mouse1/mouse1_180414';... %1
    'D:/Data/AxoAxo/mouse1/mouse1_180415';... % 2
    'D:/Data/AxoAxo/mouse1/mouse1_180501a';...%2
    'D:/Data/AxoAxo/mouse1/mouse1_180501b';...%1
    'D:/Data/AxoAxo/mouse1/mouse1_180502a';...%4
    'D:/Data/AxoAxo/mouse1/mouse1_180502b';...%0
    'D:/Data/AxoAxo/mouse3/mouse3_180627';... %0
    'D:/Data/AxoAxo/mouse3/mouse3_180628';...% 0
    'D:/Data/AxoAxo/mouse3/mouse3_180629';...%6
    'D:/Data/AxoAxo/mouse4/mouse4_181114b';...%16
    'D:/Data/AxoAxo/mouse5/mouse5_181112b';...%0
    'D:/Data/AxoAxo/mouse5/mouse5_181116';...%0
    'D:/Data/AxoAxo/mouse6/mouse6_190330';...%0
    'D:/Data/AxoAxo/mouse6/mouse6_190331';...%5
    'D:\Data\Axoaxonic_Data_Lianne\u19_200310_135409';...%3 - but should be 1 AAC, only 2? 6 and 22 definitely not
    %     'D:\Data\Axoaxonic_Data_Lianne\u19_200313_120452';...%2 % should be 4 and 41 - now 4 and 25? Pulses are not registered correctly
    'D:\Data\Axoaxonic_Data_Lianne\u19_200313_155505'};...%3 % correct - 3,13,61
    %     'D:\Data\Axoaxonic_Data_Lianne\u21_200305_153604';...
%     'D:\Data\Axoaxonic_Data_Lianne\u21_200309_142534';...
%     'D:\Data\Axoaxonic_Data_Lianne\u26_200306_172032'};...
%%

for iSess = [1:5, 17]%[16:length(dirN)]
    
    cd(dirN{iSess})
    basepath = cd; basename = bz_BasenameFromBasepath(cd);
    
    % % % % % % % % % % % % %
    % % Get AACs
    % % % % % % % % % % % % %
    [~, ~, aacs] = splitCellTypes(basepath);
    
    load([basename '.spikes.cellinfo.mat'])
    load([basename '.cell_metrics.cellinfo.mat'])
    
    load([basename '.optoStim.manipulation.mat'])
    load([basename '.dblZeta20_100_100.mat'])
    
    
    %     figure
    
    for iAAC = aacs
        cd(dirN{iSess})
        
        figure
        
        % [binnedHisto,histoEdges, spkInstPhaseStruct]=spkInstPhase(basepath, iAAC);
        [binnedHisto,histoEdges, spkInstPhaseStruct]=spkInstPhaseOutPulse(basepath, iAAC);
        
        
        subplot(2,1,1)
        
        plot(spkInstPhaseStruct.phase{1},spkInstPhaseStruct.amp{1},'.')
        xlabel('inst phase')
        ylabel('inst amp')
        %
        subplot(2,1,2)
        histogram('BinEdges', histoEdges, 'BinCounts', binnedHisto)
        hold on
        plot(histoEdges,.25*max(binnedHisto)+cos(histoEdges)*.25*max(binnedHisto),'k')
        xlabel('inst phase')
        ylabel('count')
        %     title('Preferred Theta Phase Histogram')
        
        h = suptitle({[num2str(iSess) '_' num2str(iAAC)]});
        h.Interpreter = 'none';
        
        %         savefig(gcf,['E:\Dropbox\PD_Hpc\Progress\AAC\InstAmpPhase\InstAmpPhase_AAC_' num2str(iSess) '_' num2str(iAAC) '.fig'])
        
        cd('E:\Dropbox\PD_Hpc\Progress\AAC\InstAmpPhase')
        fileStr = ['InstAmpPhase_AAC_' num2str(iSess) '_' num2str(iAAC)];
        print(gcf,[fileStr '.pdf'],'-dpdf','-bestfit')
        append_pdfs([basename '_InstAmpPhase_per_unit.pdf'],[fileStr '.pdf'])
        delete([fileStr '.pdf'])
        close
        
    end
end
