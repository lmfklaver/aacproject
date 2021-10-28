%% FR Box plots to Opto stim
% NB Split for ChR and Arch
% include lines in stead of jittered dots

launchDirNforAACSessions

%%
frAAC = []; frINT = []; frPyr = [];
frAACOpto = []; frINTOpto = []; frPyrOpto = [];
x_AAC       = []; x_INT      = []; x_Pyr       = [];
x_AACOpto    = []; x_INTOpto   = []; x_PyrOpto    = [];

for iSess = sessions
    
    basepath = dirN{iSess};
    cd(basepath)
    
    basename = bz_BasenameFromBasepath(basepath);
    
    load([basename '.spikes.cellinfo.mat']);
    load([basename '.optoStim.manipulation.mat']);
    
    pulseEpochs = optoStim.timestamps;
    timwin = [-1 1];
    binSize = 0.001;
    
    rate = getRatesTrialsBaseStim(spikes, pulseEpochs, ...
        'timwin',timwin,'binSize',binSize);
    
    baserate = rate.base;
    stimrate = rate.stim;
    
    
    meanBR = zeros(1,length(spikes.UID));
    meanSR = zeros(1,length(spikes.UID));
    for iUnit =1:length(spikes.UID)
        meanBR(iUnit) = mean(baserate{iUnit});
        meanSR(iUnit) = mean(stimrate{iUnit});
    end
    clear pyrs ints_n ints_w aacs
    
    [pyrs, ints, aacs] = splitCellTypes(basepath);
    
    frAAC = [frAAC meanBR(aacs)];
    frINT = [frINT meanBR(ints)];
    frPyr = [frPyr meanBR(pyrs)];
    
    frAACOpto = [frAACOpto meanSR(aacs)];
    frINTOpto = [frINTOpto meanSR(ints)];
    frPyrOpto = [frPyrOpto meanSR(pyrs)];
    
end

%%
% jitter scatter
x_AAC = ones(1,length(frAAC)); % cond 2
x_INT = repmat(3,1,length(frINT)); % cond 1
x_Pyr = repmat(7,1,length(frPyr));
x_AACOpto = repmat(2,1,length(frAACOpto)); % cond 2
x_INTOpto = repmat(4,1,length(frINTOpto)); % cond 1
x_PyrOpto = repmat(8,1,length(frPyrOpto));

% for jitter
x_INT  = x_INT + rand(1,length(x_INT))*0.1;
x_AAC   = x_AAC + rand(1,length(x_AAC))*0.1;
x_Pyr   = x_Pyr + rand(1,length(x_Pyr))*0.1;
x_INTOpto  = x_INTOpto + rand(1,length(x_INTOpto))*0.1;
x_AACOpto   = x_AACOpto + rand(1,length(x_AACOpto))*0.1;
x_PyrOpto   = x_PyrOpto + rand(1,length(x_PyrOpto))*0.1;

fig = figure;
set(gcf,'PaperOrientation','Landscape')

plot(x_AAC,frAAC,'o')
hold on
plot(x_INT,frINT,'o')
plot(x_Pyr,frPyr,'o')


plot(x_AACOpto,frAACOpto,'o')
hold on
plot(x_INTOpto,frINTOpto,'o')
plot(x_PyrOpto,frPyrOpto,'o')

Grps = {'frAAC' 'frAACOpto' 'frINT' 'frINTOpto'  'frPyr' 'frPyrOpto'};
GrpLbl = {'AAC','AAC Opto','INT','INT Opto','PYR','PYR Opto'};
numGrps = length(Grps);

set(gca,'XTick',1:numGrps,'XTickLabel',GrpLbl)
xlim([0.5 6.5])



allFR = [];
allFRGrp = [];

for iGrp = 1:numGrps
    allFRGrp = [allFRGrp repmat(iGrp,1,length(eval(Grps{iGrp})))];
    allFR = [allFR eval(Grps{iGrp})];
end


boxplot(allFR, allFRGrp)

set(gca,'XTick',1:numGrps,'XTickLabel',GrpLbl)
xlim([0.5 6.5])
ylabel('FR')

% title('Baseline FR')