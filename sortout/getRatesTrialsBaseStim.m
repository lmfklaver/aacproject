function [rate] = getRatesTrialsBaseStim(spikes, pulseEpochs, varargin)
% This function
%
%   USAGE 
%
%   Dependencies: 
%   Buzcode, Englishlab\utilities
%
%   INPUTS
%   Name-value paired inputs:
%   'basename'      - basename (default: basename from basepath)
%   'saveMat'       - (default: true)
%   'binSize'       -  in seconds(default: 0.01)
%   'pulsew'        - over which the stim rate is calculated, in seconds
%                     (default: 0.1)
%   'timwin'        - in seconds (default ([-0.1 0.1])
%   'stimWin'       - over how long do we want to calculate the stimrate
%                         (from 0, in seconds)
%
%   OUTPUTS
%   rate
%   .base           - rate per bin for baseline window (timwin(1) : 0)
%   .stim           - rate per bin for pulse window (0:pulsew)
%   .opts           - parameters for calculating rate
%  
%   EXAMPLES
%
%
%   NOTES
%  
%
%   TO-DO
%   - timwin(2) is redundant now, maybe also include PETH output for the
%   entire timwin
%
%   HISTORY
%   2020/12  Lianne is documenting this thing 
%
%

%% Parse! 

if ~exist('basepath','var')
    basepath = pwd;
end

basename = bz_BasenameFromBasepath(basepath);

p = inputParser;
addParameter(p,'basename',basename,@isstr);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'binSize',0.01,@isnumeric);
addParameter(p,'pulsew',0.1,@isnumeric);
addParameter(p,'timwin',[-0.1 0.1],@isnumeric);
addParameter(p,'stimWin',[ 0 0.1],@isnumeric);


parse(p,varargin{:});
basename     = p.Results.basename;
saveMat      = p.Results.saveMat;
binSize      = p.Results.binSize;
pulsew       = p.Results.pulsew;
timwin       = p.Results.timwin;
stimWin      = p.Results.stimWin;

cd(basepath)

%%

[spkTimesPerTrial] = getSpkTimTrials2(basepath, pulseEpochs);

for iUnit = 1:length(spikes.UID)
    baseInd{iUnit} = cellfun(@(a) a<0, spkTimesPerTrial{iUnit},'UniformOutput', false);
    stimInd{iUnit} = cellfun(@(a) a>0 & a <pulsew, spkTimesPerTrial{iUnit},'UniformOutput', false);
end

for iUnit = 1:length(spikes.UID)
    for iTr = 1:length(pulseEpochs)
        
        baseIndTr = baseInd{iUnit}{iTr};
        stimIndTr = stimInd{iUnit}{iTr};
        
        selUnitTr = spkTimesPerTrial{iUnit}{iTr};
        selUnitBaseTr = selUnitTr(baseIndTr);
        selUnitStimTr = selUnitTr(stimIndTr);
        
        basetrials{iUnit}{iTr}=selUnitBaseTr;
        stimtrials{iUnit}{iTr}=selUnitStimTr;
        
        baserate_temp{iUnit}{iTr} = length(basetrials{iUnit}{iTr})*1/abs(timwin(1));
        %assuming that timwin1 is negative to 0 , baseline calculated over entire timwin negative
        stimrate_temp{iUnit}{iTr} = length(stimtrials{iUnit}{iTr})*1/pulsew;%
    end
end

for iUnit = 1:length(spikes.UID)
    baserate{iUnit} = cell2mat(baserate_temp{iUnit})';
    stimrate{iUnit} = cell2mat(stimrate_temp{iUnit})';
end
   
rate.base = baserate;
rate.stim = stimrate;
rate.opts.timwin = timwin;
rate.opts.pulsew = pulsew;
rate.opts.binSize = binSize;

if saveMat
    save([basename '.rates.analysis.mat'],'rate')
end
end
    