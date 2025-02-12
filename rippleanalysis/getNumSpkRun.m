function [runspikes] = getNumSpkRun(basepath, varargin)
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
addParameter(p,'units','all',unitsValidation);

parse(p,varargin{:});
sampleRate      = p.Results.sampleRate;
saveMat         = p.Results.saveMat;
units           = p.Results.units;

cd(basepath)
%%
cd(basepath)
basename = bz_BasenameFromBasepath(basepath);
runepochs=load([basename '.run.states.mat']);
load([basename '.spikes.cellinfo.mat'],'spikes');
load([basename '.optoStim.manipulation.mat']);
pulseEpochs = optoStim.timestamps;


if strcmpi(units,'all')
    units = 1:length(spikes.times);
end

[pulseinrun] = findPulseInRun(runepochs, pulseEpochs);
[gd_run,indices] = SubtractIntervals(runepochs.run.epochs,pulseEpochs);
runtotaltime=sum(diff(gd_run'));
runpulsetotaltime=sum(diff(pulseEpochs(pulseinrun,:)'));

numSpkperRip_ONper=[];
numSpkperRip_OFFper=[];
for iUnit = units
    [status,interval] = InIntervals(spikes.times{iUnit},pulseEpochs(pulseinrun,:));
    [statusNO,intervalNO] = InIntervals(spikes.times{iUnit},gd_run); %gdrips changed from RipOutsidePulse
    RunParticipationON(iUnit)=sum(status);
    RunParticipationOFF(iUnit)=sum(statusNO);
    RunRateOn(iUnit)=sum(status)/runpulsetotaltime;
    RunRateOff(iUnit)=sum(statusNO)/runtotaltime;
end



% ripspikes.numSpkPerCycPerRip = numSpkPerCycPerRip;
runspikes.RunspikesON=RunParticipationON;
runspikes.RunspikesOFF=RunParticipationOFF;
runspikes.gd_run=gd_run;
runspikes.pulseinrun=pulseinrun;
runspikes.gd_runtime=runtotaltime;
runspikes.pulse_runtime=runpulsetotaltime;
runspikes.RunOnRate=RunRateOn;
runspikes.RunOffRate=RunRateOff
if saveMat
    save([basename '.runspikes.analysis.mat'],'runspikes');
end