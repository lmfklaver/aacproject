function ripmod = getRipMod(basepath,spikes, varargin)

%   USAGE
%   [ripmod] = getRipMod(basepath,<options>)
%
%
%   DEPENDENCIES
%
%
%
%   INPUTS
%   Name-Value pairs
%   'basename'          -
%   'saveMat'           -
%   'ccg'               -
%   'ccgBaseBins'       -
%   'binsAroundPeak'    -
%   'epochs'            -
%   'baseTime'
%
%
%   OUTPUTS
%   ripmod
%   .mod
%   .time
%
%   EXAMPLES
%  ripmod = getRipMod(basepath, spikes, 'epochs', gd_eps, 'ccg', ripple_ccg,'baseTime',[-0.4 -0.3],'baselineAroundPeak',[-.05 .05],'savemat',true)
%
%
%   NOTES
%
%
%
%   HISTORY
%   2020/12/05  Lianne pulled this out of the getPhaseMap

%   2020/12/08  Lianne documented this poorly
%   2021/11/08  Earl softcoded ripmod, and ripmod is based off of CCG
%   instead of CrossCor

%   2020/12/08  Lianne documented this

%
%   TO-DO
%
%
%

%% Parse !
if ~exist('basepath','var')
    basepath = pwd;
end

basename = bz_BasenameFromBasepath(basepath);
load([basename '.ripples.events.mat'])
p = inputParser;
addParameter(p,'basename',basename,@isstr);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'ccg',[],@isstruct);

addParameter(p,'baselineAroundPeak',[-.05 .05],@isnumeric);

addParameter(p,'epochs',[],@isnumeric);
addParameter(p, 'baseTime', [-0.5 -0.3], @isnumeric); % [start stop]

parse(p,varargin{:});
basename        = p.Results.basename;
saveMat              = p.Results.saveMat;
baselineAroundPeak  = p.Results.baselineAroundPeak;
ripple_ccg           = p.Results.ccg;
gd_eps          = p.Results.epochs;
baseTime        = p.Results.baseTime;

cd(basepath)
%%
if isempty(ripple_ccg)

    [ripple_ccg] = getRipCCG(basepath,spikes,'epochs',gd_eps,'ccgbin', 0.01,'saveMat',true); %change 10msbind

end

ripmod = [];
cellsInCCG      = 1:(size(ripple_ccg.ccg,3)-1); % if you did [{spikes} {t}]
rippleInCCG     = (size(ripple_ccg.ccg,3)); %ripple is last entry

binsAroundPeak= ripple_ccg.t >=baselineAroundPeak(1) & ripple_ccg.t <= baselineAroundPeak(2);
ccgBaseBins = ripple_ccg.t >=baseTime(1) & ripple_ccg.t <= baseTime(2);

ripmod.mod  = nanmean(ripple_ccg.ccg(binsAroundPeak,cellsInCCG,rippleInCCG),1)...
    ./nanmean(ripple_ccg.ccg(ccgBaseBins,cellsInCCG,rippleInCCG),1);
ripmod.meanbaseline=nanmean(ripple_ccg.ccg(ccgBaseBins,cellsInCCG,rippleInCCG),1);
ripmod.max = max(ripple_ccg.ccg(binsAroundPeak,cellsInCCG,rippleInCCG));
ripmod.min = min(ripple_ccg.ccg(binsAroundPeak,cellsInCCG,rippleInCCG));
ripmod.time = ripple_ccg.t;
% gd_eps=get_gd_eps(basepath);
% [status]=InIntervals(ripples.peaks(:,1),gd_eps);
% rippledur=sum(ripples.timestamps(status,2)-ripples.timestamps(status,1));
% gd_time=sum(diff(gd_eps'))-rippledur;
% for i=1:length(spikes.times);
%     [status1]=InIntervals(spikes.times{i},ripples.timestamps(status,:));
%     ripFR(i)=sum(status1)/rippledur;
%     blFR(i)=sum(~status1)/gd_time;
% end
% ripmod.modnorm=ripFR./blFR;
load([basename '.ripplepeth.analysis.mat']);
for ii=cellsInCCG
    ripmod.modPETH(ii)=nanmean(ripplepeth.count(ii,(binsAroundPeak))/(baselineAroundPeak(2)-baselineAroundPeak(1))/nanmean(ripplepeth.count(ccgBaseBins)/(baseTime(2)-baseTime(1))));
end
if saveMat
save([basename '.ripmod.mat'],'ripmod');
end
end

