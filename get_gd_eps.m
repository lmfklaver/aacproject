function gd_eps = get_gd_eps(basepath, varargin)
%
% seconds between starts of stimulations to determine stimfree epoch
%
%

if ~exist('basepath','var')
    basepath = pwd;
end

basename = bz_BasenameFromBasepath(basepath);

p = inputParser;
addParameter(p,'basename',basename,@isstring);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveAs','.gd_eps.mat',@isstring);
addParameter(p,'secToBaseline',1,@isnumeric);

parse(p,varargin{:});
basename        = p.Results.basename;
saveMat         = p.Results.saveMat;
saveAs          = p.Results.saveAs;
secToBaseline   = p.Results.secToBaseline;

cd(basepath)
%%
load([basename '.optoStim.manipulation.mat']);
lfp = bz_GetLFP('all');
st = optoStim.timestamps;

% This then finds "good episodes" not sure how good this section is, also
kp1 = find(diff([0;st(:,1)])>secToBaseline); % Finds the difference between all the start timestamps of column 1

if any(kp1>1)
    if any(kp1==1)
        gd_eps = [0 st(1,1); st(kp1(2:end)-1,2) st(kp1(2:end),1);st(end,2) inf];
    else
        gd_eps = [st(kp1-1,2) st(kp1,1);st(end,2) inf];
    end
elseif kp1==1
    gd_eps = [0 st(1,1)];
end
if gd_eps(end)==inf
   gd_eps(end)=lfp.duration
else
    if gd_eps(end)~=lfp.duration
        endofwindow=[optoStim.timestamps(end), lfp.duration]
        gd_eps=[gd_eps;endofwindow]
    end
end
if saveMat
    save([basename saveAs],'gd_eps')
end
end


