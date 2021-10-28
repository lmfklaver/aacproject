% CellExplorer
%% Subplot 
figure, 

subplot(2,1,1)

val_cond1 = cumul_ripmod_aac;
val_cond2 = cumul_ripmod_pyr;
val_cond3 = cumul_ripmod_ints;


edges = 0:0.1:3;
condition1 = histcounts(val_cond1,edges);
condition2 = histcounts(val_cond2,edges);
condition3 = histcounts(val_cond3,edges);
sumc1 = sum(condition1);
sumc2 = sum(condition2);
sumc3 = sum(condition3);


normalized1 = condition1./sumc1*100;
b1= bar(normalized1)
hold on
normalized2 = condition2./sumc2*100;
b2 = bar(normalized2)
hold on
normalized3 = condition3./sumc3*100;
b3 = bar(normalized3)

box off
xlabel('Ripple Modulation Index')
ylabel('% Occurance')
title('Ripmod Dist CellExplorer')

legend({'AAC','PYR','INTs'})
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
subplot(2,1,2)
val_cond1 = cumul_thetamod_aac;
val_cond2 = cumul_thetamod_pyr;
val_cond3 = cumul_thetamod_ints;


edges = -0.5:0.1:0.5;
condition1 = histcounts(val_cond1,edges);
condition2 = histcounts(val_cond2,edges);
condition3 = histcounts(val_cond3,edges);
sumc1 = sum(condition1);
sumc2 = sum(condition2);
sumc3 = sum(condition3);


normalized1 = condition1./sumc1*100;
b1= bar(normalized1)
hold on
normalized2 = condition2./sumc2*100;
b2 = bar(normalized2)
hold on
normalized3 = condition3./sumc3*100;
b3 = bar(normalized3)

box off
xlabel('Theta Modulation Index')
ylabel('% Occurance')
title('Thetamod Dist CellExplorer')

legend({'AAC','PYR','INTs'})
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