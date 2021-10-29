function STP = ShortTermPlasticity_axax_bare(basepath,varargin)
%
% This function is designed to get monosynaptic connections for the cells
% and find which cells are indeed significantly modulated by the
% stimulation The output of this function is used for getPhaseMap and
% Connectivity_Map
%
%
%   USAGE
%   run CellExplorer first, so there's a spikes.cluID field
%
%   % Dependencies %
%   Requires buzcode functions
%   Requires an .ait file to read optoStim from % make with makePulseFile
%
%   INPUTS
%   basepath    - path where 
%   
%   Name-value pairs:
%   'basename'  - only specify if other than basename from basepath
%   'saveMat'   - saving the results to [basename,
%                   '.burstMizuseki.analysis.mat']
%   'saveAs'    - if you want another suffix for your save
%
%   OUTPUTS
%   STP         - struct
%       .mono_con_axax_idx 
%       .mono_con_axax 
%       .pre_idx
%       .post_idx
%       .prob_uncor
%       .acg_pre
%       .acg_post
%       .prob
%       .gd_eps
%       .ses
%       .kp                 - indices of AACs
%       connection 
%
%   EXAMPLE
%
%
%   HISTORY
%   Code originially written by Sam Mckenzie
%   2020/10 Code addapted/ edited by Kaiser Arndt (7/12/20)
%   2021/02 Lianne documented and edited this code further
%   2021/05 Lianne comments
%   2021/10 Lianne replaced reading from .ait to reading from
%   optostim.manipulation.mat , and made default to find all the
%   connections
%
%   TO-DO
%   - Make kp definition and input dependent on the optogenetic manipulation
%   - Doesnt only have to be for axax only? 

%% Parse!
% 
if ~exist('basepath','var')
    basepath = pwd;
end

basename = bz_BasenameFromBasepath(basepath);

p = inputParser;
addParameter(p,'basename',basename,@isstring);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveAs','.STP.mat',@isstring);

parse(p,varargin{:});
basename        = p.Results.basename;
saveMat         = p.Results.saveMat;
saveAs          = p.Results.saveAs;

cd(basepath)

%%
gd_eps = get_gd_eps(basepath);

% this needs a cluID  
mono_res =  bz_GetMonoSynapticallyConnected(basepath,'epoch',gd_eps,'plot',false);

%% select which units to calculate shorttermplacity over

% ratemod = (nanmean(ccg(:,53:end),2) - nanmean(ccg(:,1:50),2))./nanmean(ccg(:,1:50),2);
% Question for sam on how this calculation/ why calculated this way vs how in getoptmod
%(indexing should also be softcoded)(check light artifact)

% if regexp(mod,'inhibition')
    [~,~,kp] = splitCellTypes(basepath);  % DON'T YOU JUST WANT TO CALCULATE THIS OVER ALL AND THEN SELECT?
% elseif regexp(mod,'excitation')
    
%     kp = find(optmod.p<.01 & optmod.ratemod>1);
% end

% this needs to be conditional for the type of modulation from the experiment
% >1 for excitation and <1 for inhibition . now an input, see if making a
% metadata struct and read from that - would be ideal.
% kp is the index to the cells that are axax cells that are modulated


%% define monosynaptically connected pairs  (currently being edited, not super sure what happens here)
ii=1; % start counter
mono_con_axax_idx   = []; % added LK, mono_con_axaxidx was not matching mono_res.sig_con
mono_con_axax       = [];
pre_idx             = [];
post_idx            = [];
prob_uncor          = [];
acg_pre             = [];
acg_post            = [];
prob                = []; % reserve spaces

for i = 1:length(mono_res) % uncommented, lianne
    if ~isempty(mono_res(i).sig_con)% uncommented, lianne
        mono_con_axax_idx =  ismember(mono_res.sig_con(:,2),kp); 
        % find which interneurons of the monosynaptically connected pairs are AACs,
        %  this uses cell that are modulated by stimulation, not sure if this is the best way for selecting cells

        mono_con_axax = mono_res.sig_con(mono_con_axax_idx,:); % list of AACs with monosynaptic connections
        pre_idx     = [pre_idx; mono_res.completeIndex(mono_con_axax(:,1),:)]; % find complete cell index of the presynaptic cell to the AACs
        post_idx    = [post_idx; mono_res.completeIndex(mono_con_axax(:,2),:)]; % find the complete index to the AACs
       
        for k = 1:size(mono_con_axax,1)
            prob_uncor(ii,:)    = mono_res.prob_noncor(:,mono_con_axax(k,1),mono_con_axax(k,2));
            acg_pre(ii,:)       = mono_res.prob_noncor(:,mono_con_axax(k,1),mono_con_axax(k,1)); % ACG of presynaptic cells
            acg_post(ii,:)      = mono_res.prob_noncor(:,mono_con_axax(k,2),mono_con_axax(k,2)); % ACG of postsynaptic AACs
            prob(ii,:)          = mono_res.prob(:,mono_con_axax(k,1),mono_con_axax(k,2)); % connection probably from from the mono_res
            ii= ii+1;
        end
        
        
        spikes = bz_GetSpikes('basepath', basepath);
        
        %% needed for plotting
        for i = 1:size(mono_con_axax,1)
            if size(mono_con_axax,1) == 0
                fprintf('no monosynaptic connections \n')
                continue
            else
                ref = spikes.times{mono_con_axax(i,1)};
                target = spikes.times{mono_con_axax(i,2)};
                [status] = InIntervals(ref,gd_eps);
                ref = ref(status);
                [status] = InIntervals(target,gd_eps);
                target = target(status);
                ses(i) = ShortTermCCG(ref,target,.0008,.2,'time', [ logspace(log10(5),log10(3000),20)/1000 inf]);
            end
        end
        
        %% Organize output
        if size(mono_con_axax,1) == 0
            fprintf('no monosynaptic connections \n')
            STP = 0
        else
            STP.mono_con_axax_idx = mono_con_axax_idx;
            STP.mono_con_axax = mono_con_axax;
            STP.pre_idx = pre_idx;
            STP.post_idx = post_idx;
            STP.prob_uncor = prob_uncor;
            STP.acg_pre = acg_pre;
            STP.acg_post = acg_post;
            STP.prob = prob;
            STP.gd_eps = gd_eps;
            STP.ses = ses; % really only needed for plotting

        end
    end
    STP.kp = kp;
    
    %%
    %save all variable to output file
    if saveMat
            save([basename  saveAs], 'STP')
    end
    
end

