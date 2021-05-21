% to get all the mean values for ripcycle spike count and total rip spike
% count for AACs

rcCount = 0;
ripSpkCount = 0;

for iSess = [1,2,4,5]
    
    cd(dirN{iSess})
    basepath = cd; basename = bz_BasenameFromBasepath(basepath);
    
    load([basename '.ripspikes.analysis.mat'])
    
    [~,~,aacs] = splitCellTypes(basepath);
    
    for iAAC = aacs
        if iSess == 2 && (iAAC == 10 || iAAC ==22)
            continue
        else
            rcCount = rcCount+1;
            ripSpkCount = ripSpkCount+1;
                        
            ripCycleSpk(rcCount) = mean(cell2mat(numSpkPerCycPerRip{iAAC}));
            ripSpk(ripSpkCount) = mean(spikesRipNum{iAAC});
        end
    end
end

