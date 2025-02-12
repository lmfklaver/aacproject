
function [RipRateCorr]=sm_RipRateCorr(basepath,varargin)

%   USAGE
%   Pull our ripple firing rate correlations to determine strength of
%   assemblies in and out of stimulation
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




parse(p,varargin{:});

basename    = p.Results.basename;
saveMat     = p.Results.saveMat;

load([basename '.spikes.cellinfo.mat']);
load([basename '.ripples.events.mat']);
load([basename '.optoStim.manipulation.mat']);
load([basename '.cell_metrics.cellinfo.mat']);

fils = basepath;
col = linspecer(10,'jet');
close all


max_inVin =[];
max_in_v_out =[];
max_outVout =[];
allRate =[];


ispyr =cellfun(@(a) isstr(a) && contains(a,'Pyr'),cell_metrics.putativeCellType);

    
    %%Find Firing Rate of Each Cell in Each Ripple
    for ii = 1:length(spikes.times)
        binSpk(:,ii) = cellfun(@(a) sum(spikes.times{ii}>a(1) & spikes.times{ii}<a(2))/(a(2)-a(1)),num2cell(ripples.timestamps,2)) ;       
    end
% Identify ripples during stim time
    kp = ripples.peaks>optoStim.timestamps(1) &  ripples.peaks<optoStim.timestamps(end);    
%Identify ripple peaks within optostim    
    in = InIntervals(ripples.peaks,optoStim.timestamps);
%Only spike rates in ripples during stim
    binSpk = binSpk(kp,:);
%in is all ripples during stim time, 1s being ripples within stim epochs 
    in = in(kp);

    
    
%%     
    
    %        [~,in1] = InIntervals(optoStim.timestamps(:,1),rip1(kp,:));
    %        [~,in2] = InIntervals(optoStim.timestamps(:,1),rip2(kp,:));
    %
    %        in1 = in1(in1>0);
    %        in2 = in2(in2>0);
    
    %
    %
    %          binSpk1 = binSpk1(kp,:);
    %
    %          binSpk2 = binSpk2(kp,:);
    %
    %
    %
    %         %ok = tsne(zscore(binSpk(:,ispyr)));
    %
    %         dat = zscore([binSpk1(:,ispyr);binSpk2(:,ispyr)]);
    %
    %         for ii = 1:length(in1)
    %
    %            %get corr 1st vs 2nd
    %            rr(ii) = corr(dat(in1(ii),:)',dat(in1(ii)+size(binSpk1,1),:)');
    %         end
    
    %Mean rate of all cells in each ripple, including cells without spikes
    uRate = nanmean(binSpk(:,:),2);
    [~,b] = histc(uRate,prctile(uRate,0:2:100));
    % [~,b] = histc(ripples.peaks(kp),prctile(ripples.peaks(kp),0:2:100));
    col = linspecer(50,'jet');
    %         h = figure;
    %
    % for ii = 1:50
    % plot(ok(b==ii & ~in,1),ok(b==ii & ~in,2),'.','color',col{ii},'markersize',20)
    % hold on
    % end
    %         hold on
    %
    %         plot(ok(in,1),ok(in,2),'o','color','k')
    %          plot(ok(in,1),ok(in,2),'x','color','k')
    %         print(h, '-dpsc2',filenameps ,'-append');
    %         close all
    
    
    %get nearest neighbors
    % Get only pyrs & zscore, zscored spike
    pop = zscore(binSpk(:,ispyr));
    % Correlation of pyramidal cell spike rate in ripples during opto time,
    % how correlated rates are across all ripples?
    ok = corr(pop');
    %Build grid identifying x and y as ripples
    %that occur within stim
    [XX,YY] = meshgrid(in);
    %Build template the size of xx and yy where there is no unity line in
    %center
    up = ~(eye(size(XX))==1);
    %Start calling on only correlations of interest
    
    %Correlation matrix of occurance of ripples in stim as logical?
    in_v_in =  (XX == YY) & up& (XX==1 & YY==1);   
    %In v In rate correlation
    in_v_inP = ok;
    %Make in v in rate correlation of ripples not in stim -inf to minimize
    %effect on max/mean?
    in_v_inP(~in_v_in) = -inf;
    %Only take rows/columns of ripples that happen in stim
    in_v_inP = in_v_inP(in,in);
    
    %Sort by max correlation for each ripple?
    max_inVinT =[];
    for ii = 1:size(in_v_inP,1)
        tmp = [in_v_inP(ii,:)];
        
        tmp = sort(tmp);
        max_inVinT(ii) = mean(tmp(end-9:end));
    end
    
    %concat across sessions
    max_inVin = [max_inVin;max_inVinT(:)];
    
    
    %Logical Where X or Y=1 but not both to compare in v out
    in_v_out =  (XX ~= YY) & up;
    %Entire Rate correlation
    in_v_outP = ok;
    %Select only rate correlations of interest, making others -inf
    in_v_outP(~in_v_out) = -inf;
    %Only take rows from ripples in stim
    in_v_outP = in_v_outP(in,:);
    %Sort by max correlation for each ripple?
    max_in_v_outT =[];b=[];
    for ii = 1:size(in_v_outP,1)
        tmp = [in_v_outP(ii,:)];
        %B is identified by sorting
        [tmp,bt] = sort(tmp);
        %Identify Neighbors based on max correlation
        b(ii) = bt(end);
        %Take last 10 neighbors
        max_in_v_outT(ii) = mean(tmp(end-9:end));
    end
    
    %concat across sessions
    max_in_v_out = [max_in_v_out;max_in_v_outT(:)];
    
    %all of the rows are out, how are you getting the out_ID, why not take ~in?
    out_id = b;
    out_v_out =  (XX == YY) & up & (XX==0 & YY==0);
    
    out_v_outP = ok;
    out_v_outP(~out_v_out) = -inf;
    %Sort by max correlation for each ripple?
    max_outVoutT =[];
    for ii = 1:length(out_id)
        tmp = [out_v_outP(out_id(ii),:)];
        
        [tmp] = sort(tmp);
        
        max_outVoutT(ii) = mean(tmp(end-9:end));
    end
    
    
    max_outVout = [max_outVout;max_outVoutT(:)];
    allRate = [allRate;uRate(out_id)];
   
    RipRateCorr.max_outVout=max_outVout
    RipRateCorr.allRate=allRate
    RipRateCorr.max_in_v_out=max_in_v_out
    RipRateCorr.max_inVin=max_inVin
    
    if savemat
        save([basename '.ripratecorr.analaysis.mat'])
    end
end


%ps2pdf('psfile','AAC_tsne1.ps','pdffile1','AAC_tsne.pdf','gscommand','C:\Program Files (x86)\gs\gs9.54.0\bin\gswin32.exe')
