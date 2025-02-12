function [ripple_ccgSTIM] = getRipCCGSTIM(basepath,spikes,varargin)

%   USAGE
%
%   DEPENDENCIES
%
%   INPUTS
%   basepath        -
%   spikes          -
%
%   Name-value pairs
%   'basename'      -
%   'epochs'        -
%   'saveMat'       -
%   'ccgbin'        -
%   'ccgtotsamples' -
%
%
%   OUTPUTS
%
%   EXAMPLES
%
%   NOTES
%
%   HISTORY
%
%   TO-DO
%   If no gd_eps --> gd_eps is entire session?


%% Parse !
if ~exist('basepath','var')
    basepath = pwd;
end

basename = bz_BasenameFromBasepath(basepath);

p = inputParser;
addParameter(p,'basename',basename,@isstr);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'ccgbin', 0.005,@isnumeric);
addParameter(p,'ccgtotsamples',10001,@isnumeric);
addParameter(p,'ccgdur',0.2,@isnumeric);
addParameter(p,'epochs',[],@isnumeric);



parse(p,varargin{:});

basename    = p.Results.basename;
saveMat     = p.Results.saveMat;
ccgbin      = p.Results.ccgbin;
ccgdur      = p.Results.ccgdur;
gd_eps      = p.Results.epochs;
ccgtotsamples = p.Results.ccgtotsamples;



%%

% Get ripple CCGs
cid = [];
rip_ccg = [];
NN = [];
ix = 1;

load([basename '.ripspikes.allripinstim.analysis.mat'])

rips = ripspikes.ONrips.timestamps(:,1)
selSpikes.times=[spikes.times {rips}]
[rip_ccg, t] = CCG(selSpikes.times,[],'binSize',ccgbin,'duration',ccgdur,'norm','rate');


ripple_ccgSTIM.ccg          = rip_ccg;
ripple_ccgSTIM.binsize      = ccgbin;
ripple_ccgSTIM.t            = t
ripple_ccgSTIM.ccgdur       = ccgdur;
ripple_ccgSTIM.ccglength    = ccgbin*(ccgdur); % for plotting


%%
if saveMat
    % Check if file exists:
    fripccg = [basename '.ripple_ccgSTIM.mat'];
    
    if exist(fripccg,'file')
        overwrite = input([basename,'.ripple_ccgSTIM already exists. Overwrite? [Y/N] '],'s');
        switch overwrite
            case {'y','Y'}
                delete(fripccg)
            case {'n','N'}
                return
            otherwise
                error('Y or N please...')
        end
    end
    
    save([basename '.ripple_ccgSTIM.mat'],'ripple_ccgSTIM')
end

end
