
launchDirNforAACSessions
sessions = [1:10,12];%[21 22];

allpulsepeth=[]
figure,
%% Run through sessions to collect data
% for iSess=1:length(sessions);
%     idirses=sessions(iSess);
%     cd(dirN{idirses});
%     basepath = cd;
%     basename = bz_BasenameFromBasepath(cd);
%     load ([basename '.optoStim.manipulation.mat']);
%     dirN{idirses}
%     unique(diff(optoStim.timestamps'))
% end
for iSess=1:length(sessions);
    idirses=sessions(iSess);
    cd(dirN{idirses});
    %Load session dependencies
    basepath = cd;
    basename = bz_BasenameFromBasepath(cd);
    load ([basename '_celltypes.mat']);
    load ([basename '.spikes.cellinfo.mat']);
    load ([basename '.optoStim.manipulation.mat']);
    load ([basename '.gd_eps.mat']);
    load([basename '.pulsepeth700ms10msbins.analysis.mat']);
%     [pulsepeth] = getPETH_epochs(basepath,'epochs',optoStim.timestamps(:,1),'timwin',[-.7 .7], ...
%                'binSize', 0.01,'long',true);
%     save([basename '.pulsepeth700ms10msbins.analysis.mat'], 'pulsepeth') ;
    %Threshold for PYR Cell FR below 5hz
    gd_time=sum(diff(gd_eps'));
    pyrspikes=spikes.times(pyrs);
    pyrsFR=[];
    for ipyr=1:length(pyrspikes);
        [status]=InIntervals(pyrspikes{ipyr},gd_eps);
        pyrsFR(ipyr,:)=sum(status)/gd_time;
    end
    FRThresh=find(pyrsFR'<5); %defining FR threshold
    pyrs=pyrs(FRThresh);
    allpulsepeth=[allpulsepeth; pulsepeth.rate(pyrs,:)];
end


zallpulsepeth=zscore(allpulsepeth,[],2);
sortingvals=mean(zallpulsepeth(:,80:126),2);
allpulsepethsort=[zallpulsepeth sortingvals];
sortedallpulsepeth=sortrows(allpulsepethsort,size(allpulsepethsort,2),'descend');
figure,
subplot(2,1,1)
meanpethrate=mean(zallpulsepeth,1);
bar(meanpethrate)
            ylabel('Firing Rate')
            name = {'-700';'-600';'-500';'-400';'-300';'-200';'-100';'0';'100';'200';'300';'400';'500';'600';'700'};
            title('PETH of PYR Cell FR in response to PV Inh')
set(gca,'XTick',[]);
ax = gca
ax.Position = [0.13,0.54,0.775,0.38];
ylim([-2 3]);
subplot(2,1,2);
imagesc(sortedallpulsepeth(:,1:(end-1)));
ylabel('Cell (#)');
xlabel('Time (ms)');
name = {'-700';'-600';'-500';'-400';'-300';'-200';'-100';'0';'100';'200';'300';'400';'500';'600';'700'};
xticks([0 10 20 30 40 50 60 70 80 90 100 110 120 130 140]);
set(gca,'xticklabel',name);
ax = gca;
ax.Position = [0.13,0.11,0.775,0.41];
clim([-3 4]);


