burstIndex = burstinessMizuseki_epochs(basepath,'epochs',runEpochs,'epochName', 'run5cms', 'saveMat',false);

figure,
plot(burstIndex.in, burstIndex.out,'o')
hold on
l1 = lsline;
l1.Color = 'r';
r1 = refline(1,0); % diagonal)
r1.LineStyle = '--';
xlabel('burstiness RUN')
ylabel('burstiness REST')
