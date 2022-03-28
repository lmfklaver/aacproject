function STP = ShortTermPlasticity(basepath,varargin)
%
% This function is designed to get monosynaptic connections for the cells
% and find which cells are indeed significantly modulated by the
% stimulation The output of this function is used for getPhaseMap and
% Connectivity_Map
%
%
%   USAGE
%
%
%   % Dependencies %
%   Requires buzcode functions
%   Requires the output of the function getOptoStim (ie, the struct optmod)
%
%   INPUTS
%   basepath
%
%   Name-value pairs:
%   'basename'  - only specify if other than basename from basepath
%   'saveMat'   - saving the results to [basename,
%                   '.burstMiz/useki.analysis.mat']
%   'saveAs'    - if you want another suffix for your save
%
%   OUTPUTS
%
%
%   EXAMPLE
%
%
%   HISTORY
%   Code originially written by Sam Mckenzie
%   2020/10 Code addapted/ edited by Kaiser Arndt (2020/07/12)
%   2021/02 Lianne documented and edited this code further
%   2022/03 Earl addeed code to find ACGs of non-monosynaptically connected
%   PYR in this session and monosynaptically connected in this session


%% Parse!
%
if ~exist('basepath','var')
    basepath = pwd;
end

basename = bz_BasenameFromBasepath(basepath);

kpValidation = @(x) isnumeric(x) || strcmp(x,'all');


p = inputParser;
addParameter(p,'basename',basename,@isstring);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveAs','.STP.mat',@isstring);
addParameter(p,'kp','all',kpValidation);

parse(p,varargin{:});
basename        = p.Results.basename;
saveMat         = p.Results.saveMat;
saveAs          = p.Results.saveAs;
kp = p.Results.kp;



cd(basepath)


%% Shorter version of old STP

% find gd eps
gd_eps = get_gd_eps(cd,'saveMat',false);

% Makes mono_res for only the good epochs
% what defines a good epoch
mono_res =  bz_GetMonoSynapticallyConnected(basepath,'epoch',gd_eps,'plot',false);

if strcmpi(kp,'all')
    kp = 1:size(mono_res.sig_con,1);
end

%% define monosynaptically connected pairs
ii=1; % start counter
mono_con_idx = [];
mono_con = [];
pre_idx =[];
post_idx =[];
prob_uncor =[];
acg_pre =[];
acg_post =[];
prob =[]; % reserve spaces

mono_con_idx =  ismember(mono_res.sig_con(:,2),kp);
% find which interneurons of the monosynaptically connected pairs
% are AACs, this uses cell that are modulated by stimulation, not
% sure if this is the best way for selecting cells

mono_con = mono_res.sig_con(mono_con_idx,:); %
pre_idx = [pre_idx; mono_res.completeIndex(mono_con(:,1),:)]; %
post_idx = [post_idx; mono_res.completeIndex(mono_con(:,2),:)]; %
for k = 1:size(mono_con,1)
    prob_uncor(ii,:) = mono_res.prob_noncor(:,mono_con(k,1),mono_con(k,2));
    acg_pre(ii,:)  = mono_res.prob_noncor(:,mono_con(k,1),mono_con(k,1)); % ACG of presynaptic cells
    acg_post(ii,:)  = mono_res.prob_noncor(:,mono_con(k,2),mono_con(k,2)); %
    prob(ii,:)  = mono_res.prob(:,mono_con(k,1),mono_con(k,2)); % connection probably from from the mono_res
    ii= ii+1;
end


spikes = bz_GetSpikes('basepath', basepath);
%% needed for plotting
for i = 1:size(mono_con,1)
    if size(mono_con,1) == 0
        fprintf('no monosynaptic connections \n')
        continue
    else
        
        ref = spikes.times{mono_con(i,1)};
        target = spikes.times{mono_con(i,2)};
        [status] = InIntervals(ref,gd_eps);
        ref = ref(status);
        [status] = InIntervals(target,gd_eps);
        target = target(status);
        ses(i) = ShortTermCCG(ref,target,.0008,.2,'time', [ logspace(log10(5),log10(3000),20)/1000 inf]);
    end
end

%% Organize output
if size(mono_con,1) == 0
    fprintf('no monosynaptic connections \n')
    STP.mono_con = {'no monosynaptic connections'};
else
    STP.mono_con_idx = mono_con_idx;
    STP.mono_con = mono_con;
    STP.pre_idx = pre_idx;
    STP.post_idx = post_idx;
    STP.prob_uncor = prob_uncor;
    STP.acg_pre = acg_pre;
    STP.acg_post = acg_post;
    STP.prob = prob;
    STP.gd_eps = gd_eps;
    STP.ses = ses; % really only needed for plotting
    
end
STP.kp = {kp};

%find non-monosynaptically connected PYR and monosynaptically connected Pyr
%to get average for plots
load([basename '_celltypes.mat']);
notprepyrs=[];
allpresynci=[];
notprepyrs=pyrs;
for iaac=aacs;
    presyni=find(STP.mono_con(:,2)==iaac);
    presync=STP.mono_con(presyni,1);
    [rows,cols,vals]=find(pyrs==presync);
    allpresynci=[allpresynci cols'];
end
uniquepresynci=unique(allpresynci);
notprepyrs(uniquepresynci)=[];
STP.not_mono_con_pyr=notprepyrs';
acg_notpretoaac=[]
iii=1
for iPyr=1:length(notprepyrs)
    acg_notpretoaac(iii,:)  = mono_res.prob_noncor(:,notprepyrs(iPyr),notprepyrs(iPyr)); % ACG of presynaptic cells
    iii=iii+1
end
STP.not_mono_con_to_aac_acg=acg_notpretoaac
%%
%save all variable to output file
if saveMat
    save([basename  saveAs], 'STP')
end

end
% end

