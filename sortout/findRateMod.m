function [include,p,h] = findRateMod(baseBins,stimBins,M, minsigbins)
% baseBins = baseline bins of Matrix M to calc signrank over
% stimBins = stimulation bins of Matrix M to calc signrank over
% minsigbins = minimum bins that need to be significant;
    countStim = 0;
    for iStimBin = stimBins
        countStim = countStim +1;
        countBase = 0;
        
        for iBaseBin = baseBins
            
            countBase = countBase +1;
            [p(countStim,countBase),h(countStim,countBase)] = signrank(M(:,countBase),M(:,iStimBin));
        end
    end
    include = sum(h,2)== minsigbins;
    
end

    