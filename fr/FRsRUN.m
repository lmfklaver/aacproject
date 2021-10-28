%% FR x Theta no SPW-R
% Script to FR x Theta
% addpath(genpath('E:\Dropbox\MATLAB (1)\AdrienToolBox\TStoolbox')) % Why?

launchDirNforAACSessions;

%

frAAC = []; frNINT = []; frWINT = [];   frPyr = [];
frAACRUN = []; frNINTRUN = []; frWINTRUN = []; frPyrRUN = [];
x_AAC = []; x_nINT = []; x_wINT = []; x_Pyr = [];
x_AACRUN  = []; x_nINTRUN = []; x_wINTRUN = []; x_PyrRUN = [];
sessPyr = []; sessINT = []; sessAAC = [];


for iSess = 1:length(dirN)
    
    basepath = dirN{iSess};
    cd(basepath)
    
    basename = bz_BasenameFromBasepath(basepath);
    
    % load in variables
    load([basename '.spikes.cellinfo.mat']);
    load([basename '.optoStim.manipulation.mat']);
    load([basename '.ripples.events.mat']);
    load([basename '.run.states.mat'])
    
    % 
    nonrun = zeros(1,length(run.epochs)+1)';
    nonrun(2:end,1) = run.epochs(:,2);
    nonrun(1,1) = 0;
    nonrun(1:end-1,2) = run.epochs(:,1);
    nonrun(end,:) = [];

    pulseEpochs = optoStim.timestamps;
    timwin = [-.5 .5];
    options.binSize = 0.001;
    options.stimWin = pulseEpochs(1,2)-pulseEpochs(1,1);
    
%%%NB EXCLUDE STIM PERIODS first!!! InIntervals
    [rate] = getRatesTrialsBaseStim(spikes, run.epochs);
    
    baserate= rate.base;
    runrate = rate.stim; 
%     [runrate] = getRatesTrialsBaseRUN(spikes, run.epochs);
     meanBR = []; meanSPWR =[];
    for iUnit =1:length(spikes.UID)
        meanBR(iUnit) = mean(baserate{iUnit});
        meanrun(iUnit) = mean(runrate{iUnit});
        
    end
%     [pyrs, ints_n, ints_w, aacs] = splitCellTypes_WideNarrow(basepath);
        [pyrs, ints_n,aacs] = splitCellTypes(basepath);

    length(aacs)
    frAAC   = [frAAC meanBR(aacs)];
    frNINT  = [frNINT meanBR(ints_n)];
%     frWINT  = [frWINT meanBR(ints_w)];
    frPyr   = [frPyr meanBR(pyrs)];
    
    frAACRUN    = [frAACRUN meanrun(aacs)];
    frNINTRUN   = [frNINTRUN meanrun(ints_n)];
%     frWINTRUN   = [frWINTRUN meanSPWR(ints_w)];
    frPyrRUN    = [frPyrRUN meanrun(pyrs)];
    
    sessPyr = [sessPyr repmat(iSess,1,length(frPyr))];
        sessINT = [sessINT repmat(iSess,1,length(frNINT))];
        sessAAC = [sessAAC repmat(iSess,1,length(frAAC))];
        
end
% jitter scatter
    x_AAC       = [repmat(1,1,length(frAAC))]; % cond 2
    x_nINT      = [repmat(3,1,length(frNINT))]; % cond 1
%     x_wINT      = [repmat(5,1,length(frWINT))]; % cond 3
    x_Pyr       = [repmat(5,1,length(frPyr))];
    x_AACRUN    = [repmat(2,1,length(frAAC))]; % cond 2
    x_nINTRUN   = [repmat(4,1,length(frNINT))]; % cond 1
%     x_wINTRUN   = [repmat(6,1,length(frWINT))]; % cond 3
    x_PyrRUN    = [repmat(6,1,length(frPyr))];
    
    % for jitter
    x_nINT  = x_nINT + rand(1,length(x_nINT))*0.1;
    x_AAC   = x_AAC + rand(1,length(x_AAC))*0.1;
%     x_wINT  = x_wINT + rand(1,length(x_wINT))*0.1;
    x_Pyr   = x_Pyr + rand(1,length(x_Pyr))*0.1;
    x_nINTRUN  = x_nINTRUN + rand(1,length(x_nINT))*0.1;
    x_AACRUN   = x_AACRUN + rand(1,length(x_AAC))*0.1;
%     x_wINTRUN  = x_wINTRUN + rand(1,length(x_wINT))*0.1;
    x_PyrRUN   = x_PyrRUN + rand(1,length(x_Pyr))*0.1;



%%
fig = figure;

% set(gcf,'Position',[354 634 965 704])
set(gcf,'PaperOrientation','Landscape')


plot(x_AAC,frAAC,'o')
hold on 
plot(x_nINT,frNINT,'o')
% plot(x_wINT,frWINT,'o')
plot(x_Pyr,frPyr,'o')


plot(x_AACRUN,frAACRUN,'o')
hold on 
plot(x_nINTRUN,frNINTRUN,'o')
% plot(x_wINTRUN,frWINTRUN,'o')
plot(x_PyrRUN,frPyrRUN,'o')

% Grps = {'frAAC' 'frAACRUN' 'frNINT' 'frNINTRUN' 'frWINT' 'frWINTRUN' 'frPyr' 'frPyrRUN'};
% GrpLbl = {'AAC','AAC RUN','Narrow INT','Narrow INT RUN','Wide INT','wINT RUN','PYR','PYR RUN'};
 Grps = {'frAAC' 'frAACRUN' 'frNINT' 'frNINTRUN'  'frPyr' 'frPyrRUN'};
 GrpLbl = {'AAC','AAC RUN',' INT',' INT RUN','PYR','PYR RUN'};

numGrps = length(Grps);

set(gca,'XTick',1:numGrps,'XTickLabel',GrpLbl)
xlim([0.5 6.5])

% % % % % figure, plot([frAAC ;frAACRUN],'ko-')

allFR = [];
allFRGrp = [];

for iGrp = 1:numGrps
allFRGrp = [allFRGrp repmat(iGrp,1,length(eval(Grps{iGrp})))];
allFR = [allFR eval(Grps{iGrp})];
end


boxplot([allFR], [allFRGrp])

set(gca,'XTick',1:numGrps,'XTickLabel',GrpLbl)
xlim([0.5 8.5])
ylabel('FR')


%%

 Grps = {'frAAC' 'frAACRUN' 'frNINT' 'frNINTRUN'  'frPyr' 'frPyrRUN'};
 GrpLbl = {'AAC','AAC RUN',' INT',' INT RUN','PYR','PYR RUN'};


figure, 
hold on

plot([repmat(1,1,length(frAAC));repmat(2,1,length(frAACRUN))],[frAAC ;frAACRUN],'-','Color',[211/255,211/255,211/255])
xlim([0.5 2.5])

plot([repmat(3,1,length(frNINT));repmat(4,1,length(frNINTRUN))],[frNINT ;frNINTRUN],'-','Color',[211/255,211/255,211/255])
xlim([0.5 2.5])

plot([repmat(5,1,length(frPyr));repmat(6,1,length(frPyrRUN))],[frPyr ;frPyrRUN],'-','Color',[211/255,211/255,211/255])
xlim([0.5 2.5])

allFR = [];
allFRGrp = [];

for iGrp = 1:6 %1:numGrps
allFRGrp = [allFRGrp repmat(iGrp,1,length(eval(Grps{iGrp})))];
allFR = [allFR eval(Grps{iGrp})];
end

labelsBox = {'AAC','AAC RUN','INT','INT RUN','PYR','PYR RUN'}
boxplot([allFR], [allFRGrp],'Color','km','notch','on','symbol','o','label',labelsBox)

ylabel('FR')

%%
% Histogram to match the boxplots
figure,
subplot(3,1,1)

edges = 0:5:max(frAACRUN)+10;
condition1 = histcounts(frAAC,edges);
condition2 = histcounts(frAACRUN,edges);
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
legend({'frAAC','frAACRUN'})

%
subplot(3,1,2)

edges = 0:5:max(frNINTRUN)+10;
condition1 = histcounts(frNINT,edges);
condition2 = histcounts(frNINTRUN,edges);
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
legend({'frINT','frINTRUN'})

%
subplot(3,1,3)

edges = 0:5:max(frPyrRUN)+10;
condition1 = histcounts(frPyr,edges);
condition2 = histcounts(frPyrRUN,edges);
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
legend({'frPyr','frPyrRUN'})


% % % continue