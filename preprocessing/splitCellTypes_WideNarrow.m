function [pyrs, ints_n, ints_w, aacs] = splitCellTypes_WideNarrow(basepath)

%% maybe add an session.analysisTags for Opto ChR or Arch, considering the new animals are both ChR and Arch

cd(basepath)
basename = bz_BasenameFromBasepath(cd);

load([basename '.spikes.cellinfo.mat'])
load([basename '.optoStim.manipulation.mat'])
load([basename '.cell_metrics.cellinfo.mat'])
load([basename '.dblZeta.mat'])

timwin         = [-0.5 0.5];
options.binSize = 0.001;

[baserate, stimrate] = getRatesTrialsBaseStim(spikes, optoStim.timestamps, timwin,options);

STD_Val = 3;

intsnarrow = find(contains(cell_metrics.putativeCellType,'Narrow Interneuron'));
intswide = find(contains(cell_metrics.putativeCellType,'Wide Interneuron'));

if ~isempty(regexp(basename,'mouse', 'once')) % mouse-mice are excitation
    for iUnit= 1:length(baserate)
        sd3ind(iUnit) = mean(stimrate{iUnit}) > mean(baserate{iUnit})+STD_Val*std(baserate{iUnit});
    end
    aacs = getAACnums(cell_metrics,dblZetaPChR,'excitation',sd3ind);
elseif isempty(regexp(basename,'mouse', 'once')) % u and m have been inhibition (SO FAR)
    for iUnit= 1:length(baserate)
        sd3ind(iUnit) = mean(stimrate{iUnit}) < mean(baserate{iUnit})-STD_Val*std(baserate{iUnit});
    end
    aacs = getAACnums(cell_metrics,dblZetaPArch,'inhibition',sd3ind);
end

intsindn = ~ismember(intsnarrow, aacs);
intsindw = ~ismember(intswide, aacs);
ints_n = intsnarrow(intsindn);
ints_w = intswide(intsindw);

pyrs = find(strcmpi(cell_metrics.putativeCellType,'Pyramidal Cell'));
pyrs = pyrs(~ismember(pyrs,aacs));

end
