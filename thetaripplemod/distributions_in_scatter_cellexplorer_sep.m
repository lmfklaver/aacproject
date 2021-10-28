% CellExplorer
%% Subplot 
figure, 
set(gcf,'Position',[ 1000         490         560         848])


val_cond1 = cumul_ripmod_aac;
val_cond2 = cumul_ripmod_pyr;
val_cond3 = cumul_ripmod_ints;


edges = -3:0.1:3;
condition1 = histcounts(val_cond1,edges);
condition2 = histcounts(val_cond2,edges);
condition3 = histcounts(val_cond3,edges);
sumc1 = sum(condition1);
sumc2 = sum(condition2);
sumc3 = sum(condition3);

subplot(6,1,1)

normalized1 = condition1./sumc1*100;
b1= bar(normalized1)
hold on
box off
% xlabel('Ripple Modulation Index')
ylabel('% Occurance')
title('Ripmod Dist CellExplorer')

legend({'AAC'})
forXlim = find(normalized1 ~= 0);
xlim([forXlim(1)-1 forXlim(end)+1]);
xt = get(gca,'XTick');
xl = edges(xt(2:end));
set(gca,'XTick',xt(2:end),'XTickLabel', string(xl))


subplot(6,1,2)

normalized2 = condition2./sumc2*100;
b2 = bar(normalized2);
hold on
box off
% xlabel('Ripple Modulation Index')
ylabel('% Occurance')
% title('Ripmod Dist CellExplorer')
forXlim = find(normalized2 ~= 0);
xlim([forXlim(1)-1 forXlim(end)+1]);
legend({'PYR'})
xt = get(gca,'XTick');
xl = edges(xt(2:end));
set(gca,'XTick',xt(2:end),'XTickLabel', string(xl))

subplot(6,1,3)

normalized3 = condition3./sumc3*100;
b3 = bar(normalized3)

box off
xlabel('Ripple Modulation Index')
ylabel('% Occurance')
% title('Ripmod Dist CellExplorer')
forXlim = find(normalized3 ~= 0);
xlim([forXlim(1)-1 forXlim(end)+1]);
legend({'INTs'})
xt = get(gca,'XTick');
xl = edges(xt(2:end));
set(gca,'XTick',xt(2:end),'XTickLabel', string(xl))

b1.FaceAlpha = 0.75;
b1.EdgeColor = 'none';
b1.BarWidth = 0.9;

b2.FaceAlpha = 0.75;
b2.EdgeColor = 'none';
b2.BarWidth = 0.9;

b3.FaceAlpha = 0.75;
b3.EdgeColor = 'none';
b3.BarWidth = 0.9;
%%

val_cond1 = cumul_thetamod_aac;
val_cond2 = cumul_thetamod_pyr;
val_cond3 = cumul_thetamod_ints;


edges = -3:0.1:3;
condition1 = histcounts(val_cond1,edges);
condition2 = histcounts(val_cond2,edges);
condition3 = histcounts(val_cond3,edges);
sumc1 = sum(condition1);
sumc2 = sum(condition2);
sumc3 = sum(condition3);

subplot(6,1,4)

normalized1 = condition1./sumc1*100;
b1= bar(normalized1);
hold on
box off
% xlabel('Theta Modulation Index')
ylabel('% Occurance')
title('Thetamod Dist CellExplorer')

legend({'AAC'})

forXlim = find(normalized1 ~= 0);
xlim([forXlim(1)-1 forXlim(end)+1]);

xt = get(gca,'XTick');
xl = edges(xt(2:end));
set(gca,'XTick',xt(2:end),'XTickLabel', string(xl))


subplot(6,1,5)
normalized2 = condition2./sumc2*100;
b2 = bar(normalized2);
hold on
box off
% xlabel('Theta Modulation Index')
ylabel('% Occurance')
% title('Thetamod Dist CellExplorer')

legend({'PYR'})
forXlim = find(normalized2 ~= 0);
xlim([forXlim(1)-1 forXlim(end)+1]);
xt = get(gca,'XTick');
xl = edges(xt(2:end));
set(gca,'XTick',xt(2:end),'XTickLabel', string(xl))


subplot(6,1,6)
normalized3 = condition3./sumc3*100;
b3 = bar(normalized3);

box off
xlabel('Theta Modulation Index')
ylabel('% Occurance')
% title('Thetamod Dist CellExplorer')
forXlim = find(normalized3 ~= 0);
xlim([forXlim(1)-1 forXlim(end)+1]);
legend({'INTs'})
xt = get(gca,'XTick');
xl = edges(xt(2:end));
set(gca,'XTick',xt(2:end),'XTickLabel', string(xl))

b1.FaceAlpha = 0.75;
b1.EdgeColor = 'none';
b1.BarWidth = 0.9;

b2.FaceAlpha = 0.75;
b2.EdgeColor = 'none';
b2.BarWidth = 0.9;

b3.FaceAlpha = 0.75;
b3.EdgeColor = 'none';
b3.BarWidth = 0.9;