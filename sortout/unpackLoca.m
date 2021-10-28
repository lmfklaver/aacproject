cumulCount = 0;
% sessions = [1,2,3,4,5,8,9,16,17]
getCumulRipModThetaPhase

dotSize = 50;
%%
cumulCount = 0;

for iSess = sessions
    
    cd(dirN{iSess})
    basepath    = cd;
    basename    = bz_BasenameFromBasepath(basepath);
    
    clear location
    if exist([basename '.location.mat'],'file')
        load([basename '.location.mat'])
        
        for iL = 1:length(location.iAAC)
            cumulCount = cumulCount + 1;
            cumulLocaID{cumulCount} = [num2str(iSess) '_' num2str(location.iAAC{iL})];
            cumulLocaDepth{cumulCount} = location.diffChan{iL};
        end
    else
        continue
        
    end
end


figure
% scatter(location17,cumul_ripmod_aac,dotSize,[211/255,211/255,211/255],'filled')
scatter(cell2mat(cumulLocaDepth')',abs(cumul_ripmod_aac-1),dotSize,'k','filled')

lsline
xlabel('Location')
ylabel('Ripmod')
title('Location x Ripple modulation')

[R,P] =corrcoef(abs(cumul_ripmod_aac-1),cell2mat(cumulLocaDepth'));
%P>0.05 means 