function [EpochFRCorrelation]=sm_EpochRateCorr(basepath,baselineEpoch,varargin)

%   USAGE
%   Calculate orrelation of cells that are coactive in ripples before
%   baselineEpoch, during ripples after baselineEpoch, and during ripples
%   that occur during stim after baselineEpoch, to the activity of cells
%   that are active during baselineEpoch. baselineEpoch currently
%   represents the time in a recording session that an animal is exploring
%   an environment.
%   
%  
%   DEPENDENCIES
%
%   INPUTS
%   basepath        -
%   baselineEpoch   - 2x1 time during exploration in seconds. Ripple times are removed
%   in this function
%
%   Name-value pairs
%   'basename'      -
%   'epochs'        -
%   'saveMat'       -
%
%
%   OUTPUTS
%   EpochFRCorrelation
%
%   EXAMPLES
%
%   NOTES
%
%   HISTORY
%   9/2023 Adapted from SM original code by EG


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

load([basename '.spikes.cellinfo.mat']);
load([basename '.ripples.events.mat']);
load([basename '.optoStim.manipulation.mat']);
load([basename '.cell_metrics.cellinfo.mat']);
    
i=1;
ispyr =cellfun(@(a) isstr(a) && contains(a,'Pyr'),cell_metrics.putativeCellType);
%%Find correlation of cross correlations of cells in stim ripples, pre control ripples, and post ripples?       
    CC_base_v_con_pre =[];
    CC_base_v_con_post =[];
    CC_base_v_stim =[];
    % do pairwise
    
    % get baseline spikes
    
    ep1 = optoStim.timestamps;
    %ep1(find(~(diff(ep1,[],2)==.3)),:) = [];
    ep2 =ripples.timestamps;
    ep  =MergeEpochs2([ep1;ep2]);
    %Idenify baseline relevant to the period?
    baseEp = excludeEpochs([baselineEpoch],ep);
    
    
    %find good subset of spikes
    clear spk_base
    for ii = 1:length(spikes.times)
        spk1 = spikes.times{ii};
        kp1 = InIntervals(spk1,baseEp);
        spk_base{ii} = spk1(kp1);
        
    end
    
    %find stim_rip
    kp = ripples.peaks>optoStim.timestamps(1,1) &  ripples.peaks<optoStim.timestamps(end,1);
    
    in = InIntervals(ripples.peaks,optoStim.timestamps);
    
    in = in(kp);
    gdrip = ripples.timestamps(kp,:);
    
    clear stim_rip
    for ii = 1:length(spikes.times)
        spk1 = spikes.times{ii};
        kp1 = InIntervals(spk1,gdrip(in,:));
        stim_rip{ii} = spk1(kp1);
        
    end
    
    
    clear con_rip_post
    for ii = 1:length(spikes.times)
        spk1 = spikes.times{ii};
        kp1 = InIntervals(spk1,gdrip(~in,:));
        con_rip_post{ii} = spk1(kp1);
        
    end
    
    
    pre_rip = ripples.timestamps(ripples.timestamps(:,1)<optoStim.timestamps(1,1),:);
    
    clear con_rip_pre
    for ii = 1:length(spikes.times)
        spk1 = spikes.times{ii};
        kp1 = InIntervals(spk1,pre_rip);
        con_rip_pre{ii} = spk1(kp1);
        
    end
    
    
    
    ixx=1;
    for b = logspace(log10(.0025),log10(.25),10)
        CC_base = nan(length(spikes.times));
        CC_conRip_post = nan(length(spikes.times));
        CC_conRip_pre = nan(length(spikes.times));
        CC_stimRip = nan(length(spikes.times));
        
        %CrossCorr under spike analysis - needs to be compiled
        for ii = 1:length(spikes.times)
            spk1 = con_rip_post{ii};
            if ~isempty(spk1)
                for jj = 1:length(spikes.times)
                    
                    if ii~=jj && ispyr(ii) && ispyr(jj)
                        
                        spk2 = con_rip_post{jj};
                        if ~isempty(spk2)
                            CC_conRip_post(ii,jj) = CrossCorr(spk1,spk2,b,1)/length(spk1);
                        end
                    end
                end
            end
            
        end
        
        for ii = 1:length(spikes.times)
            spk1 = con_rip_pre{ii};
            if ~isempty(spk1)
                for jj = 1:length(spikes.times)
                    
                    if ii~=jj && ispyr(ii) && ispyr(jj)
                        
                        spk2 = con_rip_pre{jj};
                        if ~isempty(spk2)
                            CC_conRip_pre(ii,jj) = CrossCorr(spk1,spk2,b,1)/length(spk1);
                        end
                    end
                end
                
                
            end
        end
        
        for ii = 1:length(spikes.times)
            spk1 = stim_rip{ii};
            if ~isempty(spk1)
                for jj = 1:length(spikes.times)
                    
                    if ii~=jj && ispyr(ii) && ispyr(jj)
                        
                        spk2 = stim_rip{jj};
                        if ~isempty(spk2)
                            
                            CC_stimRip(ii,jj) = CrossCorr(spk1,spk2,b,1)/length(spk1);
                        end
                    end
                end
                
            end
        end
        
        for ii = 1:length(spikes.times)
            spk1 = spk_base{ii};
            if ~isempty(spk1)
                for jj = 1:length(spikes.times)
                    
                    if ii~=jj && ispyr(ii) && ispyr(jj)
                        
                        spk2 = spk_base{jj};
                        if ~isempty(spk2)
                            CC_base(ii,jj) = CrossCorr(spk1,spk2,.1,1)/length(spk1);
                        end
                    end
                end
            end
            
        end
        [r,p] = corr(CC_base(:),CC_stimRip(:),'rows','pairwise');
        
        
        
        CC_base_v_stim(i,ixx) =r;
        
        [r,p] = corr(CC_base(:),CC_conRip_post(:),'rows','pairwise');
        
        
        CC_base_v_con_post(i,ixx) = [r];
        
        [r,p] = corr(CC_base(:),CC_conRip_pre(:),'rows','pairwise');
        
        
        CC_base_v_con_pre(i,ixx) = [r];
        bins(i,ixx)=b;
        
        
        ixx = ixx+1;
    end
    EpochFRCorrelation.CC_conRip_pre   =CC_conRip_pre;
    EpochFRCorrelation.CC_baseRip      =CC_base;
    EpochFRCorrelation.CC_conRip_post  =CC_conRip_post;
    EpochFRCorrelation.CC__stimRip     =CC_stimRip; 
    EpochFRCorrelation.CCbins=bins;
    EpochFRCorrelation.CC_base_v_con_pre= CC_base_v_con_pre;
    EpochFRCorrelation.CC_base_v_con_post= CC_base_v_con_post;
    EpochFRCorrelation.CC_base_v_stim= CC_base_v_stim;
    EpochFRCorrelation.CCbins=bins;
    if saveMat
    save([basename '.EpochFRCorrelation.analysis.mat'],'EpochFRCorrelation');
    end
