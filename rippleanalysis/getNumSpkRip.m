function [ripspikes] = getNumSpkRip(basepath, varargin)
%
%
%
%   USAGE
%   [spikesRipNum, numSpkPerCycPerRip] = getNumSpkRip(basepath, <options>)
%
%
%   INPUTS
%   basepath        - path that contains spikes.cellinfo.mat and
%                      ripples.events.mat
%
%   Name-value paired inputs:
%   'ripfrequency'  - Default [100 250]
%   'sampleRate'    - Default 1250
%   'saveMat'       - Default true
%   'units'         - Default 'all', optional numerical (e.g. AAC index)

%   OUTPUT
%   spikesRipNum    - number of total spikes in individual ripples
%   numSpkPerCycPerRip - number of spikes in each ripplecycle for each
%   individual ripple 
%
%   EXAMPLES
%   [spikesRipNum, numSpkPerCycPerRip] = getNumSpkRip(basepath,'units',aacs,'saveMat',true);
%
%   NOTES
%   2/2022 Added Spikes per rip code EG
%
%   TO-DO
%   - include part before first peak and after last peak
%
%%
unitsValidation = @(x) isnumeric(x) || strcmp(x,'all');


if ~exist('basepath','var')
    basepath = pwd;
end

basename = bz_BasenameFromBasepath(basepath);

p = inputParser;
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'sampleRate',1250,@isnumeric);
addParameter(p,'ripfrequency',[100 250],@isnumeric);
addParameter(p,'units','all',unitsValidation);

parse(p,varargin{:});
sampleRate      = p.Results.sampleRate;
ripfrequency    = p.Results.ripfrequency;
saveMat         = p.Results.saveMat;
units           = p.Results.units;

cd(basepath)
%%
cd(basepath)
basename = bz_BasenameFromBasepath(basepath);
load([basename '.ripples.events.mat']);
load([basename '.spikes.cellinfo.mat'],'spikes');
load([basename '.optoStim.manipulation.mat']);
load([basename '.gd_eps.mat']);
pulseEpochs = optoStim.timestamps;

if exist([basename '.ripspikes.allripinstim.analysis.mat'], 'file') == 2 % check if file exists
   load([basename '.ripspikes.allripinstim.analysis.mat']);
end

[status,interval]=InIntervals(ripples.peaks(:,1),gd_eps);
ripstart=ripples.timestamps(:,1);
ripend=ripples.timestamps(:,2);
gdrips(:,1) = ripstart(status);
gdrips(:,2) = ripend(status);
% pulls the channel from the ripples and loads the xml file
rippleChan = ripples.detectorinfo.detectionparms.channel;
lfp = bz_GetLFP(rippleChan);
lfp_ripple = BandpassFilter(double(lfp.data), sampleRate ,ripfrequency);

[~,intervalLFP] =InIntervals(lfp.timestamps, gdrips);

% find spikes in ripples
[~, interval] = cellfun(@(a) InIntervals(a, gdrips), spikes.times,'UniformOutput', false);

% clear spikesRip*
% unsure if right dimensions

spikesRipNum = cell(length(spikes.times),1);
spikesRip = cell(length(spikes.times),1);

if strcmpi(units,'all')
    units = 1:length(spikes.times);
end


% % % % Number of Spikes in the ripple

for iUnit = units
    for iRip = 1:length(gdrips)
        if ~isempty(sum(interval{iUnit}==iRip))
            spikesRip{iUnit}{iRip} = length(spikes.times{iUnit}(interval{iUnit}==iRip));
        end
    end
    spikesRipNum{iUnit} 	= cell2mat(spikesRip{iUnit});
end



    %Threshold for PYR Cell FR below 5hz
    gd_time=sum(diff(gd_eps'));

    %Need to find FR for PYR cells outside of ripples and stim for baseline
    %SWR Gain
    nonriptimes=[];
    nonrippleepoch=[];
    stimandrips=[optoStim.timestamps;ripples.timestamps]; %define intervals
    riptimesfixed=ConsolidateIntervals(stimandrips); %consolidate intervals of STIM and RIP
    nonriptimes(:,1)=[0;riptimesfixed(:,2)]; %define intervals in between previous intervals
    nonriptimes(:,2)=[riptimesfixed(:,1);lfp.duration];
    nonriptotaltime=sum(diff(nonriptimes')); %total time of good episodes
    gainbaseline=[];
    %each cells FR outside of stim and rip epochs
    for iUnit=1:length(spikes.times);
        [status]=InIntervals(spikes.times{iUnit},nonriptimes);
       gainbaseline(iUnit)=sum(status)/nonriptotaltime;
    end



% % % Number of Spikes per ripple cycle

% Now from first peak to last peak, what if spikes fall outside those peaks?
% 
% for iUnit = units
%     
%     selSpkTimes = spikes.times{iUnit};
%     cycleEp = cell(length(gdrips),1);
%     %
%     for iRip = 1:length(gdrips)
%         % pick lfp per ripple
%         selRip      = lfp_ripple(intervalLFP==iRip);
%         selRipTime  = lfp.timestamps(intervalLFP==iRip);
%         
%         % determine peaks to determine cycles
%         [~,peakInd] =  findpeaks(selRip);
%         
%         peakTime = selRipTime(peakInd);
%         numPeaks = length(peakInd);
%         
%         % make cycleEpochs
% %         clear cycleEpoch
%         cycleEpoch = cell(length(numPeaks),1);
%         
%         cycleEpoch{1} = [gdrips(iRip,1),peakTime(1)];
%         for iPk = 1:numPeaks-1
%             cycleEpoch{iPk+1}  = [peakTime(iPk) peakTime(iPk+1)];
%         end
%         
%         cycleEpoch{end+1} = [peakTime(end),gdrips(iRip,end),];
%         cycleEp{iRip}    = cell2mat(cycleEpoch);
%         
%         
%         % reshaping to get start+stops of cycle epochs
%         cycleEp{iRip}    = reshape(cycleEp{iRip},2, length(cycleEp{iRip})/2)';
%         
%         
%         % find spikes within rip cycles
%         [~, interval] = InIntervals(selSpkTimes,cycleEp{iRip});
%         
%         
%         numSpkPerCyc =[];
% %         numSpkPerCyc = cell(numPeaks-1,1);
%         
%         % see how many spikes within each cycle
%         for iPk = 1:numPeaks-1
%             numSpkPerCyc(iPk) =  length(selSpkTimes(interval==iPk));
%         end
%         
%         % and then num spks per cycle per ripple
%         numSpkPerCycPerRip{iUnit}(iRip) = sum(numSpkPerCyc); %adds all of the spikes in each cycle, doesn't differentiate
%         numSpkPerCycPerRipEach{iUnit}{iRip}=numSpkPerCyc; %shows number of spikes per rip per cycle
%     end
% end

[peakInPulse, pulseWithRip] = findRipplesInPulse(ripples, pulseEpochs);
numSpkperRip_ONper=[];
numSpkperRip_OFFper=[];
for iUnit = units
    [status,interval] = InIntervals(spikes.times{iUnit},ripples.timestamps(peakInPulse,:));
    %v = 1:length(ripples.timestamps);
    %RipOutsidePulse = ~ismember(v,RippeakInPulse);
    [statusNO,intervalNO] = InIntervals(spikes.times{iUnit},gdrips); %gdrips changed from RipOutsidePulse
    RipParticipationON{iUnit}=status;
    RipParticipationOFF{iUnit}=statusNO;
    % per ripple
    uniqueInts=unique(interval(interval~=0))';
    if isempty(uniqueInts);
    numSpkperRip_ONper(iUnit,(1:size(peakInPulse,2)))=0;
    else
    for iInterval = unique(interval(interval~=0))';
        numSpkperRip_ONper(iUnit,iInterval) = sum(length(find((interval==iInterval))));
        rateperRip_ON(iUnit,iInterval)=numSpkperRip_ONper(iUnit,iInterval)/(ripples.timestamps(iInterval,2)-ripples.timestamps(iInterval,1));
        gainperRip_ON(iUnit,iInterval)=rateperRip_ON(iUnit,iInterval)/gainbaseline(iUnit);
        end
    end
    for iIntervalNO = unique(intervalNO(intervalNO~=0))';
        numSpkperRip_OFFper(iUnit,iIntervalNO) = sum(length(find((intervalNO==iIntervalNO))));
        rateperRip_OFF(iUnit,iIntervalNO)=numSpkperRip_OFFper(iUnit,iIntervalNO)/(gdrips(iIntervalNO,2)-gdrips(iIntervalNO,1));
        gainperRip_OFF(iUnit,iIntervalNO)=rateperRip_OFF(iUnit,iIntervalNO)/gainbaseline(iUnit);
    end
    numSpkperRip_ON(iUnit) = nanmean(numSpkperRip_ONper(iUnit));
    gainRip_ON(iUnit)=(sum(numSpkperRip_ONper(iUnit,:))/sum((ripples.timestamps(peakInPulse,2)-ripples.timestamps(peakInPulse,1))))/gainbaseline(iUnit);
    if ~isempty(iIntervalNO)
        numSpkperRip_OFF(iUnit) = nanmean(numSpkperRip_OFFper(iUnit));
        gainRip_OFF(iUnit)=(sum(numSpkperRip_OFFper(iUnit,:))/sum(gdrips(:,2)-gdrips(:,1)))/gainbaseline(iUnit);
    else
        numSpkperRip_OFF(iUnit)=0;
        gainRip_OFF(iUnit)=0;
    end
        
        % gemiddelde spikes per ripple, opslaan numSpks per rip?    
end

%RippleStats
ON.peaks           = ripples.peaks(peakInPulse,:);
ON.timestamps      = ripples.timestamps(peakInPulse,:);
ON.peakNormedPower = ripples.peakNormedPower(peakInPulse,:);
ON.mean_power      = mean(ON.peakNormedPower);
ON.ripPerMin       = length(ON.peaks) / (sum((pulseEpochs(:,2)) - (pulseEpochs(:,1))) /60) ;
ON.pcDur100        = (sum(((ON.timestamps(:,2) - ON.timestamps(:,1)) * 1000)>100) / (length(ON.peaks))) * 100;
ON.pcDur80         = (sum(((ON.timestamps(:,2) - ON.timestamps(:,1)) * 1000)>80) / (length(ON.peaks))) * 100;
ON.mean_dur        = mean(ON.timestamps(:,2) - ON.timestamps(:,1))*1000;
% ON.sharpwavepeakUv = ripples.sharpwavepeakUv(peakInPulse,:);
% ON.sharpwavepeaknorm = ripples.sharpwavepeaknorm(peakInPulse,:);
% ON.sharpwavepeakZ = ripples.sharpwavepeakZ(peakInPulse,:);
% SWnumber = sum(~isnan(ripples.SW.timestamps(:,1)));
% ON.SWperc = (SWnumber / size(ON.timestamps,1));
% ON.SWpeakZScore=ripples.SW.peakZScore(peakInPulse);
% ON.SWtimestamps=ripples.SW.timestamps(peakInPulse,:);
% % collecting ripple duration in ms
% ON.rip_dur_ON = 1000*(ON.timestamps(:,2) - ON.timestamps(:,1));
% ON.SWdur=1000*(ripples.SW.timestamps(peakInPulse,2)-ripples.SW.timestamps(peakInPulse,1));

[OFFstatus,OFFinterval,OFFindex] = InIntervals(ripples.peaks(:,1),gd_eps);
logi=ones(size(ripples.peaks));
logi(peakInPulse)=0;
logi=logical(logi);
OFF.peaks           = ripples.peaks(OFFstatus);
OFF.timestamps      = ripples.timestamps(OFFstatus,:);%OFFstatus == (1),:);
OFF.peakNormedPower = ripples.peakNormedPower(OFFstatus,:);%((OFFstatus == (1)));
OFF.mean_power      = mean(OFF.peakNormedPower);
OFF.ripPerMin       = length(OFF.peaks) / (sum((gd_eps(:,2)) - (gd_eps(:,1))) /60) ;
OFF.pcDur100        = (sum(((OFF.timestamps(:,2) - OFF.timestamps(:,1)) * 1000)>100) / (length(OFF.peaks))) * 100;
OFF.pcDur80         = (sum(((OFF.timestamps(:,2) - OFF.timestamps(:,1)) * 1000)>80) / (length(OFF.peaks))) * 100;
OFF.mean_dur        = mean(OFF.timestamps(:,2) - OFF.timestamps(:,1))*1000;
% OFF.sharpwavepeakUv = ripples.sharpwavepeakUv(OFFstatus== 1);
% OFF.sharpwavepeaknorm = ripples.sharpwavepeaknorm(OFFstatus== 1);
% OFF.sharpwavepeakZ = ripples.sharpwavepeakZ(OFFstatus== 1);
% OFF.SWperc = (SWnumber / size(OFF.timestamps,1));
% OFF.SWpeakZScore=ripples.SW.peakZScore(OFFstatus,:);
% OFF.SWtimestamps=ripples.SW.timestamps(OFFstatus,:);
% % collecting ripple duration in ms
% OFF.rip_dur_OFF = 1000*(OFF.timestamps(:,2) - OFF.timestamps(:,1));
% OFF.SWdur=1000*(ripples.SW.timestamps(OFFstatus,2)-ripples.SW.timestamps(OFFstatus,1));


%ripspikes.numSpkPerCycPerRipEach=numSpkPerCycPerRipEach;
%ripspikes.numSpkPerCycPerRip = numSpkPerCycPerRip;
ripspikes.spikesRipNum = spikesRipNum;
ripspikes.numSpkperRip_ON=numSpkperRip_ONper;
ripspikes.numSpkperRip_OFF=numSpkperRip_OFFper;
ripspikes.avgSpkPerRip_ON=numSpkperRip_ON;
ripspikes.avgSpkPerRip_OFF=numSpkperRip_OFF;
ripspikes.RipParticipation_ON= RipParticipationON;
ripspikes.RipParticipation_OFF= RipParticipationOFF;
ripspikes.rateperRip_ON=rateperRip_ON;
ripspikes.gainperRip_ON=gainperRip_ON;
ripspikes.gainRip_ON=gainRip_ON;
ripspikes.rateperRip_OFF=rateperRip_OFF;
ripspikes.gainperRip_OFF=gainperRip_OFF;
ripspikes.gainRip_OFF=gainRip_OFF;
ripspikes.ONrips=ON
ripspikes.OFFrips=OFF
if saveMat
    save([basename '.ripspikes.allripinstim.analysis.mat'],'ripspikes');
end
