

basepath = cd; basename = bz_BasenameFromBasepath(cd);

load([basename '.ripples.events.mat'])
load([basename '.spikes.cellinfo.mat'])

[status, interval] = cellfun(@(a) InIntervals(a, ripples.timestamps), spikes.times,'UniformOutput', false);
% spikesRipTest = cellfun(@(a,b) a(b), spikes.times, status);

for iUnit = 1:length(spikes.times)
    
    for iRip = 1:length(ripples.timestamps)
        if ~isempty(sum(interval{iUnit}==iRip))
%             spikesRip{iUnit}{iRip} = length(spikes.times{iUnit}(interval{iUnit}==iRip));
            spikesRip{iUnit}{iRip} = (spikes.times{iUnit}(interval{iUnit}==iRip));
        end
    end
    
    spikesRipNum{iUnit} 	= cell2mat(spikesRip{iUnit});
    rippleSpkRate{iUnit}    = spikesRipNum{iUnit}'./(ripples.timestamps(:,2)-ripples.timestamps(:,1));
end


[pyrs, ints, aacs] = splitCellTypes(basepath);

edgesRip = 0:11;

for iAAC = aacs
    figure
    
    histogram(spikesRipNum{iAAC},edgesRip)
    xlabel('Num Spk in Rip')
    ylabel('count')
    set(gca,'YScale','log','XTick',[0.5:1:10.5],'XTickLabel', num2cell(edgesRip))
end

%% convert this to rate, also rate outside ripple
fils  = getAllExtFiles(basepath,'rip',1);
rip = LoadEvents(fils{1});

% pulls the channel from the ripples and loads the xml file
rippleChan = str2num(rip.description{1}(regexp(rip.description{1},'[0-9]')));


lfp = bz_GetLFP(rippleChan);
lfp_ripple = BandpassFilter(double(lfp.data), 1250 ,[100 250]);

%%

%%

[~,intervalLFP] =InIntervals(lfp.timestamps, ripples.timestamps);

for iAAC = aacs
    selSpkTimes = spikes.times{iAAC};
    %%
    for iRip = 1:length(ripples.peaks)
        
        % pick lfp per ripple
        selRip      = lfp_ripple(intervalLFP==iRip);
        selRipTime  = lfp.timestamps(intervalLFP==iRip);
        
        % determine peaks to determine cycles
        [~,peakInd] =  findpeaks(selRip);
        
        peakTime = selRipTime(peakInd);
        numPeaks = length(peakInd);
        
        % make cycleEpochs
        clear cycleEpoch
        cycleEpoch{1} = [ripples.timestamps(iRip,1),peakTime(1)];
        for iPk = 1:numPeaks-1
            cycleEpoch{iPk+1}  = [peakTime(iPk) peakTime(iPk+1)];
        end
        
        cycleEpoch{end+1} = [peakTime(end),ripples.timestamps(iRip,end),];
        
        cycleEp{iRip}    = cell2mat(cycleEpoch);
        
        
        % reshaping to get start+stops of cycle epochs
        cycleEp{iRip}    = reshape(cycleEp{iRip},2, length(cycleEp{iRip})/2)';
        
        
        
        % find spikes within rip cycles
        [status, interval] = InIntervals(selSpkTimes,cycleEp{iRip});
        
        
        clear numSpkPerCyc
        
        % see how many spikes within each cycle
        for iPk = 1:numPeaks-1
            numSpkPerCyc(iPk) =  length(selSpkTimes(interval==iPk));
        end
        
        % and then num spks per cycle per ripple
        numSpkPerCycPerRip{iRip} = numSpkPerCyc;
    end
            
                edgesRip = 0:1:5;
        figure,
        histogram(cell2mat(numSpkPerCycPerRip),edgesRip)
        xlabel('Avg Number of Spikes per Ripple Cycle')
        ylabel('count')
        set(gca,'YScale','log','XTick',edgesRip,'XTickLabel', num2cell(edgesRip))
        box off
%         clear numSpkPerCycPerRip
        
end


