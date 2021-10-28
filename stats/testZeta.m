
addpath(genpath('C:\Users\lklaver\Documents\GitHub\ZETA\'))

baseWinBefore = 0.3;
totalDurWin = 0.6;

for iUnit = 1:length(spikes.UID), 
    [dblZetaP(iUnit),vecLatencies(iUnit,:), sZETA{iUnit},sRate{iUnit}] = getZeta(spikes.times{iUnit},pulseEpochs-baseWinBefore,totalDurWin,[],4);
end