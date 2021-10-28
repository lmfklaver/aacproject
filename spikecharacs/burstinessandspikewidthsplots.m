% scatters

sessions = 1,2,4,5,15,16,17;
%%
pullBurstinessIndxfromSpikes

%%
depth = [12.5
50
25
0
175
-50
0
-25
25
0
-62.35
0
-137.5
-112.5
0
-137.5
-112.4];

% 1,2,4,5,15,16,17;
%%

fig =figure,
subplot(1,3,1)

% depth = a(:,6)
scatter(depth,a(:,3),'filled')
title('Depth x Burstiness,Royer')

subplot(1,3,2)
scatter(depth,a(:,4),'filled')
title('Depth x Burstiness,Mizuseki')

subplot(1,3,3)
scatter(depth,a(:,5),'filled')
title('Depth x Burstiness,Doublets')

han = axes(fig,'visible','off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han, 'Burstiness Index')
xlabel(han, 'Distance from center PYR layer')
%%

fig =figure,
% subplot(1,3,1)

% depth = a(:,6)
scatter(depth,spkwidth(:,3),'filled')
title('Depth x SpikeWidth')

ylabel('SpikeWidth')
xlabel('Distance from center PYR layer')

%%
getCumulRipModThetaPhase  % needs sessions as an input

%%
fig =figure,
% subplot(1,3,1)

% depth = a(:,6)
scatter(depth,cumul_ripmod_aac,'filled')
title('Depth x Ripmod')

ylabel('Ripmod')
xlabel('Distance from center PYR layer')
cumul_ripmod_aac