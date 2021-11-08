function [ripple_ccg] = getRipCCG(basepath,spikes,varargin)

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
%      [ripple_ccg] = getRipCCG(basepath,spikes,'epochs',gd_eps,'ccgbin', 0.001,'saveMat',true);
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
addParameter(p,'ccgbin', 0.001,@isnumeric);
addParameter(p,'ccgdur',1,@isnumeric);
addParameter(p,'epochs',[],@isnumeric);



parse(p,varargin{:});
basename    = p.Results.basename;
saveMat     = p.Results.saveMat;
ccgbin      = p.Results.ccgbin;
ccgdur      = p.Results.ccgdur;
gd_eps      = p.Results.epochs;



%%
if isempty(gd_eps)
    disp('No gd_eps')
    % OR: gd_eps is entire session
end

% Get ripple CCGs
cid = [];
rip_ccg = [];
NN = [];
ix = 1;

load([basename '.ripples.events.mat'])
load([basename '.spikes.cellinfo.mat'])

t = ripples.peaks; %11/5/2021 changed ripple onset to ripple peaks

for j = 1:length(spikes.times)
    [status] = InIntervals(spikes.times{j},gd_eps);
    
    if isfield(spikes,'cluID')
        cid = [cid;spikes.shankID(j) spikes.cluID(j)];
    else
        cid = [cid;spikes.shankID(j) spikes.UID(j)];
    end
    
    spikesGd{j} = (spikes.times{j}(status));
end
spikesandrip = [spikesGd,{t}];
    %CCG takes cell arrays of Nx1
    [rip_ccg,ccgtime] = CCG(spikesandrip,[],'Fs',30000,'binSize',ccgbin,...
        'duration',ccgdur,'norm','rate');
    
   

ripple_ccg.ccg          = rip_ccg;
ripple_ccg.t             =ccgtime;
ripple_ccg.binsize      = ccgbin;
ripple_ccg.ccgdur       = ccgdur;
ripple_ccg.ccglength     = ccgbin*(ccgdur); % for plotting

%%
if saveMat
    % Check if file exists:
    fripccg = [basename '.ripple_ccg.mat'];
    
    if exist(fripccg,'file')
        overwrite = input([basename,'.ripple_ccg already exists. Overwrite? [Y/N] '],'s');
        switch overwrite
            case {'y','Y'}
                delete(fripccg)
            case {'n','N'}
                return
            otherwise
                error('Y or N please...')
        end
    end
    
    save([basename '.ripple_ccg.mat'],'ripple_ccg')
end

end
