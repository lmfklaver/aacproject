function [pyrs, ints, aacs] = splitCellTypes(basepath)
% This function
%
%   USAGE 
%
%   Dependencies: 
%   Buzcode, Englishlab\utilities,spikes.cellinfo.mat,
%   optoStim.manipulation.mat, cell_metrics.cellinfo.mat,
%   .dblZeta20_100_100.mat
%
%   INPUTS

%   Name-value paired inputs:
%
%   OUTPUTS
%
%
%   EXAMPLES
%
%
%   NOTES
%  
%
%   TO-DO
%   - maybe add an session.analysisTags for Opto ChR or Arch, considering the new animals are both ChR and Arch
%   - Move stats + getRatesTrialsBaseStim into the AAC functions   
%
%   HISTORY
%% 

cd(basepath)
basename = bz_BasenameFromBasepath(cd);

load([basename '.spikes.cellinfo.mat'])
load([basename '.optoStim.manipulation.mat'])
load([basename '.cell_metrics.cellinfo.mat'])

% load([basename '.dblZeta.mat'])
% if no ZETA, calculate ZETA? 
load([basename '.dblZeta20_100_100.mat'])

[baserate, stimrate] = getRatesTrialsBaseStim(spikes, optoStim.timestamps, 'timwin', [-0.5 0.5],'binSize', 0.001,'pulsew',0.05);

STD_Val = 2;

if ~isempty(regexp(basename,'mouse', 'once')) % mouse-mice are excitation
    dblZetaPChR     = dblZetaPChR100;
    
    for iUnit= 1:length(baserate)
        InclusionInd(iUnit) = mean(stimrate{iUnit}) > mean(baserate{iUnit});%+STD_Val*std(baserate{iUnit});
    end
    intsall = find(contains(cell_metrics.putativeCellType,'Interneuron'));
    aacs = getAACnums(cell_metrics,dblZetaPChR,'excitation',InclusionInd);
    intsind = ~ismember(intsall, aacs);
    ints = intsall(intsind);
    
elseif isempty(regexp(basename,'mouse', 'once')) % u and m have been inhibition (SO FAR)
    dblZetaPArch    = dblZetaPArch100;
    
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
