
%% Subplot 
figure, 
set(gcf,'Position',[ 1000         490         560         848])

val_cond1 = cumul_ripmod_aac;
val_cond2 = cumul_ripmod_pyr;
val_cond3 = cumul_ripmod_ints;


edges = 0:0.1:2;
condition1 = histcounts(val_cond1,edges);
condition2 = histcounts(val_cond2,edges);
condition3 = histcounts(val_cond3,edges);
sumc1 = sum(condition1);
sumc2 = sum(condition2);
sumc3 = sum(condition3);

subplot(6,1,1)
normalized1 = condition1./sumc1*100;
h1= bar(normalized1);
hold on
box off
% xlabel('Ripple Modulation Index')
ylabel('% Occurance')
title('Ripmod Dist Samcode')

legend({'AAC'},'Location','northwest')
xt = get(gca,'XTick');
xl = edges(xt(2:end));
set(gca,'XTick',xt(2:end),'XTickLabel', string(xl))


subplot(6,1,2)
normalized2 = condition2./sumc2*100;
h2 = bar(normalized2);
box off
% xlabel('Ripple Modulation Index')
ylabel('% Occurance')
% title('Ripmod Dist Samcode')
legend({'PYR'},'Location','northwest')
xt = get(gca,'XTick');
xl = edges(xt(2:end));
set(gca,'XTick',xt(2:end),'XTickLabel', string(xl))

subplot(6,1,3)
normalized3 = condition3./sumc3*100;
h3 = bar(normalized3);
box off
xlabel('Ripple Modulation Index')
ylabel('% Occurance')
% title('Ripmod Dist Samcode')

legend({'INTs'},'Location','northwest')
xt = get(gca,'XTick');
xl = edges(xt(2:end));
set(gca,'XTick',xt(2:end),'XTickLabel', string(xl))

h1.FaceAlpha = 0.75;
h1.EdgeColor = 'none';
h1.BarWidth = 0.9;

h2.FaceAlpha = 0.75;
h2.EdgeColor = 'none';
h2.BarWidth = 0.9;

h3.FaceAlpha = 0.75;
h3.EdgeColor = 'none';
h3.BarWidth = 0.9;
%%
val_cond1 = cumul_thetaphase_aac;
val_cond2 = cumul_thetaphase_pyr;
val_cond3 = cumul_thetaphase_ints;


edges = -pi:0.1:pi;
condition1 = histcounts(val_cond1,edges);
condition2 = histcounts(val_cond2,edges);
condition3 = histcounts(val_cond3,edges);
sumc1 = sum(condition1);
sumc2 = sum(condition2);
sumc3 = sum(condition3);

subplot(6,1,4)
normalized1 = condition1./sumc1*100;
% h1= histogram(normalized1,edges);
b1 = bar(normalized1)
hold on
box off
% xlabel('Theta Phase')
ylabel('% Occurance')
title('Thetaphase Dist Samcode')
% xlim([-pi pi])

legend({'AAC'},'Location','northwest')

%
subplot(6,1,5)
normalized2 = condition2./sumc2*100;
% h2 = histogram(normalized2,edges);
b2 = bar(normalized2)
box off
% xlabel('Theta Phase')
ylabel('% Occurance')
% title('Thetaphase Dist Samcode')

legend({'PYR'},'Location','northwest')

% xlim([-pi pi])
% hold on

%
subplot(6,1,6)
normalized3 = condition3./sumc3*100;
% h3 = histogram(normalized3,edges)
b3 = bar(normalized3)

box off
xlabel('Theta Phase')
ylabel('% Occurance')
% title('Thetaphase Dist Samcode')

legend({'INTs'},'Location','northwest')
% xlim([-pi pi])

b1.FaceAlpha = 0.75;
b1.EdgeColor = 'none';
b1.BarWidth = 0.9;

b2.FaceAlpha = 0.75;
b2.EdgeColor = 'none';
b2.BarWidth = 0.9;

b3.FaceAlpha = 0.75;
b3.EdgeColor = 'none';
b3.BarWidth = 0.9;