%% FR x SPW-R

% To do:
% - ripples outside pulses, outside ripples
% - check rate function


launchDirNforAACSessions;

%%
frAAC = []; frINT = []; frPyr = [];
frAACRip = []; frINTRip = []; frPyrRip = [];
sessPyr = [];       sessINT = [];    sessAAC = [];
x_AAC       = []; x_nINT      = [];  x_Pyr       = [];
x_AACRip    = [];    x_nINTRip   = []; x_PyrRip    = [];


for iSess = sessions
    
    basepath = dirN{iSess};
    cd(basepath)
    basename = bz_BasenameFromBasepath(basepath);
    
    % load the variables
    load([basename '.spikes.cellinfo.mat']);
    load([basename '.optoStim.manipulation.mat']);
    load([basename '.ripples.events.mat']);
    
    pulseEpochs     = optoStim.timestamps;
    timwin          = [-.5 .5];
    binSize         = 0.001;
    stimWin         = pulseEpochs(1,2)-pulseEpochs(1,1);
    
    rate        = getRatesTrialsBaseStim(spikes, pulseEpochs, 'timwin', timwin,'binSize',binSize,'stimWin',stimWin);
    baserate    = rate.base;
    riprate     = getRatesTrialsBaseRip(spikes, ripples);
    
    meanBR = zeros(1,length(spikes.UID));
    meanSPWR = zeros(1,length(spikes.UID));
    
    for iUnit =1:length(spikes.UID)
        meanBR(iUnit)   = mean(baserate{iUnit});
        meanSPWR(iUnit) = mean(riprate{iUnit});
    end
    
    [pyrs, ints,aacs] = splitCellTypes(basepath);
    
    frAAC   = [frAAC meanBR(aacs)];
    frINT  = [frINT meanBR(ints)];
    frPyr   = [frPyr meanBR(pyrs)];
    
    frAACRip    = [frAACRip meanSPWR(aacs)];
    frINTRip   = [frINTRip meanSPWR(ints)];
    frPyrRip    = [frPyrRip meanSPWR(pyrs)];
    
    sessPyr = [sessPyr repmat(iSess,1,length(frPyr))];
    sessINT = [sessINT repmat(iSess,1,length(frINT))];
    sessAAC = [sessAAC repmat(iSess,1,length(frAAC))];
    
end


% jitter scatter
x_AAC       = [ones(1,length(frAAC))]; % cond 2
x_INT      = [repmat(3,1,length(frINT))]; % cond 1
x_Pyr       = [repmat(5,1,length(frPyr))];
x_AACRip    = [repmat(2,1,length(frAAC))]; % cond 2
x_INTRip   = [repmat(4,1,length(frINT))]; % cond 1
x_PyrRip    = [repmat(6,1,length(frPyr))];

% for jitter
x_INT  = x_INT + rand(1,length(x_INT))*0.1;
x_AAC   = x_AAC + rand(1,length(x_AAC))*0.1;
x_Pyr   = x_Pyr + rand(1,length(x_Pyr))*0.1;
x_INTRip  = x_INTRip + rand(1,length(x_INT))*0.1;
x_AACRip   = x_AACRip + rand(1,length(x_AAC))*0.1;
x_PyrRip   = x_PyrRip + rand(1,length(x_Pyr))*0.1;



%% And plot
% Boxplot jittered dots and boxes

fig = figure;

% set(gcf,'Position',[354 634 965 704])
set(gcf,'PaperOrientation','Landscape')


plot(x_AAC,frAAC,'o')
hold on
plot(x_INT,frINT,'o')
plot(x_Pyr,frPyr,'o')


plot(x_AACRip,frAACRip,'o')
hold on
plot(x_INTRip,frINTRip,'o')
plot(x_PyrRip,frPyrRip,'o')

Grps = {'frAAC' 'frAACRip' 'frINT' 'frINTRip'  'frPyr' 'frPyrRip'};
GrpLbl = {'AAC','AAC Rip',' INT',' INT Rip','PYR','PYR Rip'};

numGrps = length(Grps);

set(gca,'XTick',1:numGrps,'XTickLabel',GrpLbl)
xlim([0.5 6.5])

% % % % % figure, plot([frAAC ;frAACRip],'ko-')

allFR = [];
allFRGrp = [];

for iGrp = 1:numGrps
    allFRGrp = [allFRGrp repmat(iGrp,1,length(eval(Grps{iGrp})))];
    allFR = [allFR eval(Grps{iGrp})];
end


boxplot(allFR, allFRGrp)

set(gca,'XTick',1:numGrps,'XTickLabel',GrpLbl)
xlim([0.5 8.5])
ylabel('FR')


%%
% Boxplots lines and boxes


Grps = {'frAAC' 'frAACRip' 'frINT' 'frINTRip'  'frPyr' 'frPyrRip'};
GrpLbl = {'AAC','AAC Rip',' INT',' INT Rip','PYR','PYR Rip'};


figure,
hold on

plot([ones(1,length(frAAC));repmat(2,1,length(frAACRip))],[frAAC ;frAACRip],'-','Color',[211/255,211/255,211/255])
xlim([0.5 2.5])

plot([repmat(3,1,length(frINT));repmat(4,1,length(frINTRip))],[frINT ;frINTRip],'-','Color',[211/255,211/255,211/255])
xlim([0.5 2.5])

plot([repmat(5,1,length(frPyr));repmat(6,1,length(frPyrRip))],[frPyr ;frPyrRip],'-','Color',[211/255,211/255,211/255])
xlim([0.5 2.5])

allFR = [];
allFRGrp = [];

for iGrp = 1:length(Grps) %1:numGrps
    allFRGrp = [allFRGrp repmat(iGrp,1,length(eval(Grps{iGrp})))];
    allFR = [allFR eval(Grps{iGrp})];
end

labelsBox = {'AAC','AAC Rip','INT','INT Rip','PYR','PYR Rip'};
boxplot(allFR, allFRGrp,'Color','km','notch','on','symbol','o','label',labelsBox)

ylabel('FR')

%%
% Histogram to match the boxplots
figure,
subplot(3,1,1)

edges = 0:5:max(frAACRip)+10;
condition1 = histcounts(frAAC,edges);
condition2 = histcounts(frAACRip,edges);
sumc1 = sum(condition1);
sumc2 = sum(condition2);

normalized1 = condition1./sumc1*100;
h1 = histogram('BinEdges', edges, 'BinCounts', normalized1,...
    'FaceColor','k','EdgeColor','none');
h1.FaceAlpha = 0.5;
hold on
normalized2 = condition2./sumc2*100;
h2  = histogram('BinEdges', edges, 'BinCounts', normalized2,...
    'FaceColor','m','EdgeColor','none');
h2.FaceAlpha = 0.5 ;

hold on, box off
xlabel('FR'), ylabel('% Occurance'), title('AAC')
legend({'frAAC','frAACRip'})

%
subplot(3,1,2)

edges = 0:5:max(frINTRip)+10;
condition1 = histcounts(frINT,edges);
condition2 = histcounts(frINTRip,edges);
sumc1 = sum(condition1);
sumc2 = sum(condition2);

normalized1 = condition1./sumc1*100;
h1 = histogram('BinEdges', edges, 'BinCounts', normalized1,...
    'FaceColor','k','EdgeColor','none');
h1.FaceAlpha = 0.5;
hold on
normalized2 = condition2./sumc2*100;
h2  = histogram('BinEdges', edges, 'BinCounts', normalized2,...
    'FaceColor','m','EdgeColor','none');
h2.FaceAlpha = 0.5 ;

hold on, box off
xlabel('FR'), ylabel('% Occurance'), title('INT')
legend({'frINT','frINTRip'})

%
subplot(3,1,3)

edges = 0:5:max(frPyrRip)+10;
condition1 = histcounts(frPyr,edges);
condition2 = histcounts(frPyrRip,edges);
sumc1 = sum(condition1);
sumc2 = sum(condition2);

normalized1 = condition1./sumc1*100;
h1 = histogram('BinEdges', edges, 'BinCounts', normalized1,...
    'FaceColor','k','EdgeColor','none');
h1.FaceAlpha = 0.5;
hold on
normalized2 = condition2./sumc2*100;
h2  = histogram('BinEdges', edges, 'BinCounts', normalized2,...
    'FaceColor','m','EdgeColor','none');
h2.FaceAlpha = 0.5 ;

hold on, box off
xlabel('FR'), ylabel('% Occurance'), title('PYR')
legend({'frPyr','frPyrRip'})

%%
 Grps = {'frAAC' 'frAACRip' 'frNINT' 'frNINTRip'  'frPyr' 'frPyrRip'};
 GrpLbl = {'AAC','AAC Rip',' INT',' INT Rip','PYR','PYR Rip'};


figure,
hold on

plot([ones(1,length(frAAC));repmat(2,1,length(frAACRip))],[frAAC ;frAACRip],'-','Color',[211/255,211/255,211/255])
xlim([0.5 2.5])

plot([repmat(3,1,length(frNINT));repmat(4,1,length(frNINTRip))],[frNINT ;frNINTRip],'-','Color',[211/255,211/255,211/255])
xlim([0.5 2.5])

plot([repmat(5,1,length(frPyr));repmat(6,1,length(frPyrRip))],[frPyr ;frPyrRip],'-','Color',[211/255,211/255,211/255])
xlim([0.5 2.5])

allFR = [];
allFRGrp = [];

for iGrp = 1:length(Grps) %1:numGrps
allFRGrp = [allFRGrp repmat(iGrp,1,length(eval(Grps{iGrp})))];
allFR = [allFR eval(Grps{iGrp})];
end

labelsBox = {'AAC','AAC Rip','INT','INT Rip','PYR','PYR Rip'};
boxplot(allFR, allFRGrp,'Color','km','notch','on','symbol','o','label',labelsBox)

ylabel('FR')

