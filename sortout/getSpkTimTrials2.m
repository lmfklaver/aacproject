function [spkTimesPerTrial] = getSpkTimTrials2(basepath, epochs, varargin)
% This function
%
%   USAGE
%
%   Dependencies:
%
%
%   INPUTS
%   
%
%   Name-value paired inputs:
%  'timwin'
%  'binSize'
%
%   OUTPUTS
%
%
%   EXAMPLES
%
%
%   NOTES
%  
%
%   TO-DO
% 
%
%   HISTORY
% 
%
%
%
%% Parse!
if ~exist('basepath','var')
    basepath = pwd;
end

basename = bz_BasenameFromBasepath(basepath);

p = inputParser;
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'timwin',[-0.4 0.4],@isnumeric);
addParameter(p,'binSize',0.01,@isnumeric);


parse(p,varargin{:});
saveMat       = p.Results.saveMat;
timwin        = p.Results.timwin;
binSize       = p.Results.binSize;

cd(basepath)
%%
load([basename '.spikes.cellinfo.mat'])
%%

timeBefore = abs(timwin(1));
timeAfter = timwin(2);

trlCenteredPulseStart = epochs(:,1)-timeBefore;
trlCenteredPulseStop = epochs(:,1)+timeAfter;

trlCenteredPulse = [trlCenteredPulseStart trlCenteredPulseStop];

%%
spike_toPulse = realignSpikes(spikes, trlCenteredPulse);

%%
for iUnit = 1:length(spikes.UID)
    plotSpkOffset = 0;
    
    for iEp = 1:length(epochs)
        spikeTrl_Pulse{iUnit}{iEp} = spike_toPulse{iUnit}{iEp} - epochs(iEp,1);
    end
end

spkTimesPerTrial = spikeTrl_Pulse;


end
    