function ripmod = getRipMod_CCG(basepath, varargin)

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
%
%
%   OUTPUTS
%   ripmod
%   .mod
%   .time
%
%   EXAMPLES
%
%
%
%   NOTES
%
%
%
%   HISTORY
%   2020/12/05  Lianne pulled this out of the getPhaseMap
%   2020/12/08  Lianne documented this
%
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

p = inputParser;
addParameter(p,'basename',basename,@isstr);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'ccg',[],@isstruct);
addParameter(p,'ccgBaseBins', 1:80,@isnumeric);
addParameter(p,'binsAroundPeak',20,@isnumeric);
addParameter(p,'epochs',[],@isnumeric);


parse(p,varargin{:});
basename        = p.Results.basename;
saveMat         = p.Results.saveMat;
ccgBaseBins     = p.Results.ccgBaseBins;
binsAroundPeak  = p.Results.binsAroundPeak;
ripple_ccg      = p.Results.ccg;
gd_eps          = p.Results.epochs;

cd(basepath)
%%
gd_eps = get_gd_eps(basepath);
if isempty(gd_eps)
    disp('No gd_eps, skipping session')
    return
else 

if isempty(ripple_ccg)
    load([basename '.spikes.cellinfo.mat'])
    [ripple_ccg] = getRipCCG_CCG(basepath,spikes,'epochs',gd_eps,'ccgbin', 0.005,'ccgtotsamples',10001,'saveMat',false);
end

t = ripple_ccg.t;

for iUnit = 1:length(spikes.times)
    ccg(:,iUnit) = ripple_ccg.ccg(:,end,iUnit);
end

ripmod = [];
ccgPeakBin  =   (length(t)-1)/2;
ripmod.mod  =    nanmean(ccg(ccgPeakBin-binsAroundPeak:ccgPeakBin+binsAroundPeak,:))./nanmean(ccg(ccgBaseBins,:));
ripmod.ripccg = ripple_ccg;
ripmod.basename  = basename;

save([basename '.ripmod.analysis.mat'],'ripmod')
end
end

