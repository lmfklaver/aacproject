function [spkTimesPerTrial] = getSpkTimTrials2(basepath, epochs, varargin)
%   This function is calculating the spike times per trial, and checks (in
%   getRatesTrialsStim), if the number of spikes in each trial before onset
%   of epochs (e.g. pulseEpochs), is larger or smaller than during the
%   epochs
%
%   USAGE
%   [spkTimesPerTrial] = getSpkTimTrials2(basepath, epochs, varargin)
%
%   Dependencies:
%   buzcode toolbox, 
%   realignSpikes
%
%   INPUTS
%   basepath    - where the spikes.cellinfo struct is located
%   epochs      - what do you want to center your spiketimes to
%
%   Name-value paired inputs:
%  'timwin' - (Default [-0.4 0.4] time before and after epoch%  onset
%
%   OUTPUTS
%   spkTimesPerTrial    - cellArray for each unit, for each trial, spike
%                           times aligned to epoch start (at 0 s)
%
%   EXAMPLES
%   [spkTimesPerTrial] = getSpkTimTrials2(basepath, pulseEpochs);
%
%   NOTES 
%
%   TO-DO
%
%   HISTORY
%       2021/05/25 LK commented + tried with stripped path. 
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


parse(p,varargin{:});
saveMat       = p.Results.saveMat;
timwin        = p.Results.timwin;

cd(basepath)
%%
load([basename '.spikes.cellinfo.mat'])
%% 
timeBefore  = abs(timwin(1));
timeAfter   = timwin(2);

epochStart  = epochs(:,1)-timeBefore;
epochStop   = epochs(:,1)+timeAfter;

trlCenteredEpoch = [epochStart epochStop];

%%
spike_toPulse = realignSpikes(spikes, trlCenteredEpoch);

%%
for iUnit = 1:length(spikes.UID)
    plotSpkOffset = 0;
    
    for iEp = 1:length(epochs)
        spkTimesPerTrial{iUnit}{iEp} = spike_toPulse{iUnit}{iEp} - epochs(iEp,1);
    end
end

end
    