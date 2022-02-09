function [spikesRipNum, numSpkPerCycPerRip] = getNumSpkRip(basepath, varargin)
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
%
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
load([basename '.spikes.cellinfo.mat'],'spikes')


fils  = getAllExtFiles(basepath,'rip',1);
rip = LoadEvents(fils{1});

% pulls the channel from the ripples and loads the xml file
rippleChan = str2double(rip.description{1}(regexp(rip.description{1},'[0-9]')));
lfp = bz_GetLFP(rippleChan);
lfp_ripple = BandpassFilter(double(lfp.data), sampleRate ,ripfrequency);

[~,intervalLFP] =InIntervals(lfp.timestamps, ripples.timestamps);

% find spikes in ripples
[~, interval] = cellfun(@(a) InIntervals(a, ripples.timestamps), spikes.times,'UniformOutput', false);

% clear spikesRip*
% unsure if right dimensions

spikesRipNum = cell(length(spikes.times),1);
spikesRip = cell(length(spikes.times),1);

if strcmpi(units,'all')
    units = 1:length(spikes.times);
end


% % % % Number of Spikes in the ripple

for iUnit = units
    for iRip = 1:length(ripples.timestamps)
        if ~isempty(sum(interval{iUnit}==iRip))
            spikesRip{iUnit}{iRip} = length(spikes.times{iUnit}(interval{iUnit}==iRip));
        end
    end
    spikesRipNum{iUnit} 	= cell2mat(spikesRip{iUnit});
end




% % % % Number of Spikes per ripple cycle

% % Now from first peak to last peak, what if spikes fall outside those peaks?

for iUnit = units
    
    selSpkTimes = spikes.times{iUnit};
    cycleEp = cell(length(ripples.timestamps),1);
    %
    for iRip = 1:length(ripples.peaks)
        
        % pick lfp per ripple
        selRip      = lfp_ripple(intervalLFP==iRip);
        selRipTime  = lfp.timestamps(intervalLFP==iRip);
        
        % determine peaks to determine cycles
        [~,peakInd] =  findpeaks(selRip);
        
        peakTime = selRipTime(peakInd);
        numPeaks = length(peakInd);
        
        % make cycleEpochs
%         clear cycleEpoch
        cycleEpoch = cell(length(numPeaks),1);
        
        cycleEpoch{1} = [ripples.timestamps(iRip,1),peakTime(1)];
        for iPk = 1:numPeaks-1
            cycleEpoch{iPk+1}  = [peakTime(iPk) peakTime(iPk+1)];
        end
        
        cycleEpoch{end+1} = [peakTime(end),ripples.timestamps(iRip,end),];
        cycleEp{iRip}    = cell2mat(cycleEpoch);
        
        
        % reshaping to get start+stops of cycle epochs
        cycleEp{iRip}    = reshape(cycleEp{iRip},2, length(cycleEp{iRip})/2)';
        
        
        % find spikes within rip cycles
        [~, interval] = InIntervals(selSpkTimes,cycleEp{iRip});
        
        
        clear numSpkPerCyc
%         numSpkPerCyc = cell(numPeaks-1,1);
        
        % see how many spikes within each cycle
        for iPk = 1:numPeaks-1
            numSpkPerCyc(iPk) =  length(selSpkTimes(interval==iPk));
        end
        
        % and then num spks per cycle per ripple
        numSpkPerCycPerRip{iUnit}(iRip) = sum(numSpkPerCyc);
    end
end

if saveMat
    save([basename '.ripspikes.analysis.mat'],'numSpkPerCycPerRip','spikesRipNum');
end