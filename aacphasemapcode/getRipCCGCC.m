function [ripple_ccg_cc] = getRipCCGCC(basepath,spikes,varargin)

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
%      [ripple_ccg_cc] = getRipCCGCC(basepath,spikes,'epochs',gd_eps,'ccgbin', 0.005,'ccgtotsamples',10001,'saveMat',true);
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
ccgtotsamples = p.Results.ccgtotsamples;
ccgdur      = p.Results.ccgdur;
gd_eps      = p.Results.epochs;



%%
if isempty(gd_eps)
    disp('No gd_eps')
    % OR: gd_eps is entire session
end


% Get ripple CCGs
cid = [];
rip_ccg_cc = [];
NN = [];
ix = 1;

load([basename '.ripples.events.mat'])
load([basename '.spikes.cellinfo.mat'])
% rip = LoadEvents([basepath '/' basename '.evt.rip']);
% t   = rip.time(cellfun(@any,regexp(rip.description,'start')));
t = ripples.timestamps(:,1);

for j = 1:length(spikes.times)
    [status] = InIntervals(spikes.times{j},gd_eps);
    
    if isfield(spikes,'cluID')
        cid = [cid;spikes.shankID(j) spikes.cluID(j)];
    else
        cid = [cid;spikes.shankID(j) spikes.UID(j)];
    end
    
    newVar{j} = (spikes.times{j}(status));
    rip_ccg_cc(ix,:) = CrossCorr(t,spikes.times{j}(status),ccgbin,ccgtotsamples)/length(t);

    
    NN(ix) = sum(status);
    ix = ix+1;
end    
    


ripple_ccg_cc.ccg          = rip_ccg_cc;
ripple_ccg_cc.binsize      = ccgbin;
ripple_ccg_cc.binlength    = ccgtotsamples;
% ripple_ccg.ccgdur       = ccgdur;
ripple_ccg_cc.ccglength     = ccgbin*(ccgtotsamples-1); % for plotting

%%
if saveMat
    % Check if file exists:
    fripccg = [basename '.ripple_ccg_cc.mat'];
    
    if exist(fripccg,'file')
        overwrite = input([basename,'.ripple_ccg_cc already exists. Overwrite? [Y/N] '],'s');
        switch overwrite
            case {'y','Y'}
                delete(fripccg)
            case {'n','N'}
                return
            otherwise
                error('Y or N please...')
        end
    end
    
    save([basename '.ripple_ccg_cc.mat'],'ripple_ccg_cc')
end

end
