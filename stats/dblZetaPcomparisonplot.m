% dblZeta P Scatter
load([basename '.dblZeta20_300.mat'])
[pyrs, ints, aacs] = splitCellTypes(basepath);


%%
legendAACs = cell(1,length(aacs)+1);
legendAACs{1,1} = 'Others';

for iAAC=1:length(aacs)
legendAACs{1,iAAC+1} = num2str(aacs(iAAC));
end



dotSize = 100;
figure


scatter(dblZetaPArch20,dblZetaPArch100,dotSize,'k')
hold on 

for iAAC=1:length(aacs)
scatter(dblZetaPArch20(aacs(iAAC)),dblZetaPArch100(aacs(iAAC)),dotSize,'filled')
end

xlabel('Zeta 20ms');
ylabel('Zeta 100ms');

title({[num2str(basename)] 'dblZeta over first 20ms vs 100 ms'})

xlim([0 1])
ylim([0 1])
axis square

xline(0.05,':')
yline(0.05,':')
legend(legendAACs)


%%
load([basename '.dblZeta20_300.mat'])
[pyrs, ints, aacs] = splitCellTypes(basepath);


%%
legendAACs = cell(1,length(aacs)+1);
legendAACs{1,1} = 'Others';

for iAAC=1:length(aacs)
legendAACs{1,iAAC+1} = num2str(aacs(iAAC));
end



dotSize = 100;
figure


scatter(dblZetaPArch20,dblZetaPArch100,dotSize,'k')
hold on 

for iAAC=1:length(aacs)
scatter(dblZetaPArch20(aacs(iAAC)),dblZetaPArch100(aacs(iAAC)),dotSize,'filled')
end

xlabel('dblZetaPArch20ms');
ylabel('dblZetaPArch100ms');

title({[num2str(basename)] 'dblZeta over first 20ms vs 300 ms'})

xlim([0 1])
ylim([0 1])
axis square

xline(0.05,':')
yline(0.05,':')
legend(legendAACs)