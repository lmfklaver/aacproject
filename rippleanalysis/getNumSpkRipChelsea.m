function [ripspikes] = getNumSpkRipChelsea(basepath, varargin)
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
%   11/2023 edited to remove optostim base for chelsea - EG
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
gdrips=ripples.timestamps;

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
gd_eps=[0 lfp.duration];

% % % % Number of Spikes in the ripple

for iUnit = units
    for iRip = 1:length(gdrips)
        if ~isempty(sum(interval{iUnit}==iRip))
            spikesRip{iUnit}{iRip} = length(spikes.times{iUnit}(interval{iUnit}==iRip));
        end
    end
    spikesRipNum{iUnit} 	= cell2mat(spikesRip{iUnit});
end

    nonriptotaltime=sum(diff([0;lfp.duration])); %total time of good episodes
    gainbaseline=[];
    %each cells FR outside of stim and rip epochs
    for iUnit=1:length(spikes.times);
       [status]=InIntervals(spikes.times{iUnit},gd_eps);
       gainbaseline(iUnit)=sum(status)/(spikes.times{iUnit}(end)-spikes.times{iUnit}(1));
    end



% % % Number of Spikes per ripple cycle

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

numSpkperRip_OFFper=[];
for iUnit = units
    [statusNO,intervalNO] = InIntervals(spikes.times{iUnit},gdrips); %gdrips changed from RipOutsidePulse
    RipParticipationOFF{iUnit}=statusNO;
    % per ripple
    for iIntervalNO = unique(intervalNO(intervalNO~=0))';
        numSpkperRip_OFFper(iUnit,iIntervalNO) = sum(length(find((intervalNO==iIntervalNO))));
        rateperRip_OFF(iUnit,iIntervalNO)=numSpkperRip_OFFper(iUnit,iIntervalNO)/(gdrips(iIntervalNO,2)-gdrips(iIntervalNO,1));
        gainperRip_OFF(iUnit,iIntervalNO)=rateperRip_OFF(iUnit,iIntervalNO)/gainbaseline(iUnit);
    end
    if ~isempty(iIntervalNO)
        numSpkperRip_OFF(iUnit) = nanmean(numSpkperRip_OFFper(iUnit));
        gainRip_OFF(iUnit)=(sum(numSpkperRip_OFFper(iUnit,:))/sum(gdrips(:,2)-gdrips(:,1)))/gainbaseline(iUnit);
    else
        numSpkperRip_OFF(iUnit)=0;
        gainRip_OFF(iUnit)=0;
    end
        
        % gemiddelde spikes per ripple, opslaan numSpks per rip?    
end

[OFFstatus,OFFinterval,OFFindex] = InIntervals(ripples.peaks(:,1),gd_eps);
OFF.peaks           = ripples.peaks(OFFstatus);
OFF.timestamps      = ripples.timestamps(OFFstatus,:);%OFFstatus == (1),:);
OFF.peakNormedPower = ripples.peakNormedPower(OFFstatus,:);%((OFFstatus == (1)));
OFF.mean_power      = mean(OFF.peakNormedPower);
OFF.ripPerMin       = length(OFF.peaks) / (sum((gd_eps(:,2)) - (gd_eps(:,1))) /60) ;
OFF.pcDur100        = (sum(((OFF.timestamps(:,2) - OFF.timestamps(:,1)) * 1000)>100) / (length(OFF.peaks))) * 100;
OFF.pcDur80         = (sum(((OFF.timestamps(:,2) - OFF.timestamps(:,1)) * 1000)>80) / (length(OFF.peaks))) * 100;
OFF.mean_dur        = mean(OFF.timestamps(:,2) - OFF.timestamps(:,1))*1000;
OFF.sharpwavepeakUv = ripples.sharpwavepeakUv(OFFstatus== 1);
OFF.sharpwavepeaknorm = ripples.sharpwavepeaknorm(OFFstatus== 1);
OFF.sharpwavepeakZ = ripples.sharpwavepeakZ(OFFstatus== 1);
OFF.SWperc = (SWnumber / size(OFF.timestamps,1));
OFF.SWpeakZScore=ripples.SW.peakZScore(OFFstatus,:);
OFF.SWtimestamps=ripples.SW.timestamps(OFFstatus,:);
% collecting ripple duration in ms
OFF.rip_dur_OFF = 1000*(OFF.timestamps(:,2) - OFF.timestamps(:,1));
OFF.SWdur=1000*(ripples.SW.timestamps(OFFstatus,2)-ripples.SW.timestamps(OFFstatus,1));


%ripspikes.numSpkPerCycPerRipEach=numSpkPerCycPerRipEach;
%ripspikes.numSpkPerCycPerRip = numSpkPerCycPerRip;
ripspikes.spikesRipNum = spikesRipNum;
ripspikes.numSpkperRip_OFF=numSpkperRip_OFFper;
ripspikes.avgSpkPerRip_OFF=numSpkperRip_OFF;
ripspikes.RipParticipation_OFF= RipParticipationOFF;
ripspikes.rateperRip_OFF=rateperRip_OFF;
ripspikes.gainperRip_OFF=gainperRip_OFF;
ripspikes.gainRip_OFF=gainRip_OFF;
ripspikes.OFFrips=OFF
% if saveMat
%     save([basename '.ripspikes.analysis.mat'],'ripspikes');
% end
