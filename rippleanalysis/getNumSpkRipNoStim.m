function [ripspikesNoStim] = getNumSpkRipNoStim(basepath, varargin)
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


numSpkperRip_OFFper=[];
for iUnit = units
    %v = 1:length(ripples.timestamps);
    %RipOutsidePulse = ~ismember(v,RippeakInPulse);
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
        
end

%%Cell explorer gain - find spike rate within 150ms surrounding ripple peak,
%%and divide by spike rate 75ms before
[status,interval]=InIntervals(ripples.peaks(:,1),gd_eps);
CERip(:,1)=ripples.peaks(status)-.150
CERip(:,2)=ripples.peaks(status)+.150
CEBL(:,1)=CERip(:,1)-.75
CEBL(:,2)=CERip(:,1)
for iUnit = units
    [~,intervalRIP] = InIntervals(spikes.times{iUnit},CERip); %gdrips changed from RipOutsidePulse
    [~,intervalBL] = InIntervals(spikes.times{iUnit},CEBL); %gdrips changed from RipOutsidePulse
    % per ripple
    for iIntervalNO = unique(intervalRIP(intervalRIP~=0))';
        numSpkperRip_OFFRIP(iUnit,iIntervalNO) = sum(length(find((intervalRIP==iIntervalNO))))/.150;
        numSpkperRip_OFFBL(iUnit,iIntervalNO) = sum(length(find((intervalBL==iIntervalNO))))/.75; %% What if no spikes in baseline? %%INF!
        gainRipCE_OFF(iUnit,iIntervalNO)=(numSpkperRip_OFFRIP(iUnit,iIntervalNO)-numSpkperRip_OFFBL(iUnit,iIntervalNO))/(numSpkperRip_OFFRIP(iUnit,iIntervalNO)+numSpkperRip_OFFBL(iUnit,iIntervalNO));
    end
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
% OFF.sharpwavepeakUv = ripples.sharpwavepeakUv(OFFstatus== 1);
% OFF.sharpwavepeaknorm = ripples.sharpwavepeaknorm(OFFstatus== 1);
% OFF.sharpwavepeakZ = ripples.sharpwavepeakZ(OFFstatus== 1);
% % SWnumber = sum(~isnan(ripples.SW.timestamps(:,1)));
% OFF.SWperc = (SWnumber / size(OFF.timestamps,1));
% OFF.SWpeakZScore=ripples.SW.peakZScore(OFFstatus,:);
% OFF.SWtimestamps=ripples.SW.timestamps(OFFstatus,:);
% collecting ripple duration in ms
OFF.rip_dur_OFF = 1000*(OFF.timestamps(:,2) - OFF.timestamps(:,1));
% OFF.SWdur=1000*(ripples.SW.timestamps(OFFstatus,2)-ripples.SW.timestamps(OFFstatus,1));


%ripspikes.numSpkPerCycPerRipEach=numSpkPerCycPerRipEach;
%ripspikes.numSpkPerCycPerRip = numSpkPerCycPerRip;
ripspikesNoStim.spikesRipNum = spikesRipNum;
ripspikesNoStim.numSpkperRip_OFF=numSpkperRip_OFFper;
ripspikesNoStim.avgSpkPerRip_OFF=numSpkperRip_OFF;
ripspikesNoStim.RipParticipation_OFF= RipParticipationOFF;
ripspikesNoStim.rateperRip_OFF=rateperRip_OFF;
ripspikesNoStim.gainperRip_OFF=gainperRip_OFF;
ripspikesNoStim.gainRip_OFF=gainRip_OFF;
ripspikesNoStim.OFFrips=OFF;
ripspikesNoStim.gainRipCE_OFF=gainRipCE_OFF;
if saveMat
    save([basename '.ripspikes.NoStim.analysis.mat'],'ripspikesNoStim');
end
