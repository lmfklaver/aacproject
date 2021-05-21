function [pyrs, ints, aacs] = splitCellTypes(basepath)
% This function splits PYRs, INTs and AACs using cellexplorer and
% optotagging as classifiers 
%
%   USAGE 
%
%   Dependencies: 
%   Buzcode, Englishlab\utilities,spikes.cellinfo.mat,
%   optoStim.manipulation.mat, cell_metrics.cellinfo.mat,
%    zeta mat file.
%
%   INPUTS
%   basepath
%   Name-value paired inputs:
%
%   OUTPUTS
%   pyrs    - indices for Pyramidal Cells
%   ints    - indices for Interneurons without AACs
%   aacs    - indices for AACs
%
%   EXAMPLES
%
%
%   NOTES
%  
%
%   TO-DO
%   - maybe add an session.analysisTags for Opto ChR or Arch, considering the new animals are both ChR and Arch
%   - allow stats input
%
%   HISTORY
 %  Updating this in AAC_PRoject, not in Utilties NB
%% 

cd(basepath)
basename = bz_BasenameFromBasepath(cd);

load([basename '.spikes.cellinfo.mat'])
load([basename '.optoStim.manipulation.mat'])
load([basename '.cell_metrics.cellinfo.mat'])

% load([basename '.dblZeta.mat'])
% if no ZETA, calculate ZETA? 
load([basename '.pethzeta.stats.mat'])


% timwin          = [-0.5 0.5];
% options.binSize = 0.001;
% options.stimWin = 0.05;
 
[rate] = getRatesTrialsBaseStim(spikes, optoStim.timestamps, 'timwin', [-0.5 0.5],'binSize', 0.001,'pulsew',0.05);
baserate = rate.base;
stimrate = rate.stim;

STD_Val = 2;

if ~isempty(regexp(basename,'mouse', 'once')) % mouse-mice are excitation
    dblZetaPChR     = zeta.P;%dblZetaPChR100;
    
    for iUnit= 1:length(baserate)
        InclusionInd(iUnit) = mean(stimrate{iUnit}) > mean(baserate{iUnit});%+STD_Val*std(baserate{iUnit});
    end
    intsall = find(contains(cell_metrics.putativeCellType,'Interneuron'));
    aacs = getAACnums(cell_metrics,dblZetaPChR,'excitation',InclusionInd);
    intsind = ~ismember(intsall, aacs);
    ints = intsall(intsind);
    
elseif isempty(regexp(basename,'mouse', 'once')) % u and m have been inhibition (SO FAR)
    dblZetaPArch    = zeta.P;%dblZetaPArch100;
    
    for iUnit= 1:length(baserate)
        InclusionInd(iUnit) = mean(stimrate{iUnit}) < mean(baserate{iUnit});%;-STD_Val*std(baserate{iUnit});
    end
    
    intsall = find(contains(cell_metrics.putativeCellType,'Interneuron'));
    aacs    = getAACnums(cell_metrics,dblZetaPArch,'inhibition',InclusionInd);
    
    intsind = ~ismember(intsall, aacs);
    ints    = intsall(intsind);
end

pyrs = find(strcmpi(cell_metrics.putativeCellType,'Pyramidal Cell'));
pyrs = pyrs(~ismember(pyrs,aacs));

end
