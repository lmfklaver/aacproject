% Below this is plotting
%%
a = ismember(STP.mono_con_axax(:,1), STP.kp)


%everything below is plotting


%%%%%%%%%%%%%%%%%

% plot stim PETH


figure
imagesc(sortby(nanzscore(ccg,[],2),sign(optmod.ratemod).*(1./p))) % what is this p? LK)
colormap('jet')

% ccg and ratemod not of equal length


%% Works


% plot all phase portraits of axax
close all
figure
ix=1;
ph_bin = linspace(-pi,pi,16);
k = gaussian2Dfilter([10 10],[.5 .5]);

kp = find(optmod.p<.01 & optmod.ratemod>1);
ax = tight_subplot(6,7);
for i = kp' % find back this kp
    axes(ax(ix))
    imagesc(ph_bin,[],nanconvn((ph_mod.ph_rate(:,1:end-1,i)),k),[min(linearize(ph_mod.ph_rate(:,1:end-1,i)))...
        max(linearize(ph_mod.ph_rate(:,1:end-1,i)))])
    hold on
    colormap('jet')
    plot(ph_mod.ph_bin,10+cos(ph_mod.ph_bin)*10,'w')
    set(gca,'ytick',0:10:length(ph_mod.freq),'yticklabel',round(ph_mod.freq(1:10:end)),'xticklabel',[])
    set(gca,'ydir','normal')

%     title(num2str(cid((i),:)))
    ix = ix+1;
end
%% Error
% plot presyn phase portraits
ph_bin = linspace(-pi,pi,16);
% kp = find(p<.01 & ratemod>1);
figure
ix=1;
k = gaussian2Dfilter([10 10],[.05 .05]);
ax = tight_subplot(6,6);
for i = 1:length(kp)
    axes(ax(ix))

    pre_kp = pre_idx(ismember(post_idx(:,[1 4]) ,axax(i,:),'rows'),:);


    imagesc(ph_bin,[],nanconvn(nanmean(ph_rate1(:,:,ismember(cid,pre_kp(:,[1 4]),'rows')),3),k))
    hold on
    colormap('jet')
    plot(ph_bin,10+cos(ph_bin)*10,'w')
    set(gca,'ytick',1:2:length(freq),'yticklabel',round(freq(1:2:end)),'xticklabel',[])
    set(gca,'ydir','normal')

    title(num2str(cid((kp(i)),:)))
    ix = ix+1;
end

%% Works


%plot spike transmission phase portraits
close all
trans_phase = cellArrayTo3D({ph_map(1:end).trans_phase});
ph_bin = linspace(-pi,pi,32);
kp = find(p<.01 & ratemod>1);
figure
ix=1;
k = gaussian2Dfilter([20 10],[2 1]);
ax = tight_subplot(6,6);
for i = 1:length(kp)
    axes(ax(ix))

    kp_con = (ismember(mono_con(:,[1 3]) ,axax(i,:),'rows'));


    imagesc(ph_bin,[],nanconvn(nanmedian(nanzscore(trans_phase(:,:,kp_con),[],2),3),k))
    %imagesc(ph_bin,[],nanconvn(nanmedian((trans_phase(:,:,kp_con)),3),k))
    hold on
    colormap('jet')
    plot(ph_bin,10+cos(ph_bin)*10,'w')
    set(gca,'ytick',1:10:length(freq),'yticklabel',round(freq(1:10:end)),'xticklabel',[])
    set(gca,'ydir','normal')

    title(num2str(cid((kp(i)),:)))
    ix = ix+1;
end


%% Works

%plot axo axo stim PETHs
figure


ax = tight_subplot(6,6);
ccg = cell2mat(cellfun(@(a) a.ccg,optmod(2:end),'uni',0)');

p = cell2mat(cellfun(@(a) a.p,optmod(2:end),'uni',0)');
oMod = cell2mat(cellfun(@(a) a.optoMod,optmod(2:end),'uni',0)');
stimrate = cell2mat(cellfun(@(a) a.stimrate,optmod(2:end),'uni',0)');
ix=1;

kp = find(p<.01 & ratemod>1);
for i = kp'
    axes(ax(ix))
    bar(ccg(i,:))
    ix = ix+1;
end

%% Works

%plot axo axo ripple PETHs
figure
ax = tight_subplot(6,6);
ccg = cell2mat(cellfun(@(a) a.ccg,optmod(2:end),'uni',0)');

p = cell2mat(cellfun(@(a) a.p,optmod(2:end),'uni',0)');
oMod = cell2mat(cellfun(@(a) a.optoMod,optmod(2:end),'uni',0)');
stimrate = cell2mat(cellfun(@(a) a.stimrate,optmod(2:end),'uni',0)');
ratemod = (nanmean(ccg(:,53:end),2) - nanmean(ccg(:,1:50),2))./nanmean(ccg(:,1:50),2);
kp = find(p<.01 & ratemod>1);
ix=1;
for i = 1:length(kp')
    axes(ax(ix))
    if idx(i) == 1 % added loop to plot diff cell types into diff colors
        plot(-25:.005:25,rip_ccg(kp(i),:),'k')
    else
        plot(-25:.005:25,rip_ccg(kp(i),:),'r')
    end
    xlim([-1 1])
    ix = ix+1;

end


set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 0.5, ...
'DefaultAxesXColor', 'k', ...
'DefaultAxesYColor', 'k', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 12, ...
'DefaultAxesFontName', 'Arial', ...
'DefaultLineLineWidth', 1, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 12, ...
'DefaultTextFontName', 'Arial', ...
'DefaultAxesBox', 'off', ...
'DefaultAxesTickLength', [0.02 0.025]);

% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');



% To restore a property to its original MATLAB® default, use the 'remove' keyword.
% set(groot,'DefaultFigureColormap','remove')

figure
plot(-25:.005:25,sum(rip_ccg(kp(idx==1),:),1)/length(kp(idx==1)),'k')
hold on


plot(-25:.005:25,sum(rip_ccg(kp(idx==2),:),1)/length(kp(idx==2)),'r')
xlim([-1 1])

box off



%% Works



% plot CCG & ACG feature for two axax types (must have important short term
% plasticity)

%plot

figure
plot(nanmean(prob_uncor(kp1,:)),'k')
hold on
plot(nanmean(prob_uncor(kp2,:)),'r')

figure
plot(nanmean(acg_pre(kp1,:)),'k')
hold on
plot(nanmean(acg_pre(kp2,:)),'r')


figure
plot(nanmean(acg_post(kp1,:)),'k')
hold on
plot(nanmean(acg_post(kp2,:)),'r')


%% Works
% plot axo axo clusters

figure
plot3(ph_pref_theta(kp(idx==1)),ph_pref_gam(kp(idx==1)),ripmod(kp(idx==1)),'o','color','k')
hold on

plot3(ph_pref_theta(kp(idx==2)),ph_pref_gam(kp(idx==2)),ripmod(kp(idx==2)),'o','color','r')

xlabel('Theta phase preference')
ylabel('Gamma phase preference')
zlabel('Ripple modulation')

%% Works

% plot short term plasticity
close all
k = gaussian2Dfilter([1 100],[.001 1]);
n = cell2mat(arrayfun(@(a) a.time.n(1,:),ses,'uni',0)');
ok = cellfun(@(a) nanPad(a,20),{ses.trans},'uni',0);

ok = cell2mat(ok');
okt =[];
for i = 1:size(ok,1)

    okt(i,:) = nanconvn(ok(i,1:end),k);
end


% plot STP curves
figure
semilogx(ses(1).time.bins(1:end-1),nanmean((okt(kp1,:))),'k')
hold on
semilogx(ses(1).time.bins(1:end-1),nanmean((okt(kp2,:))),'r')


% plot distribution of resonance peaks
figure
[~,b1] = max(okt(kp1,:),[],2);
[~,b2] = max(okt(kp2,:),[],2);
semilogx(ses(1).time.bins,histc(ses(1).time.bins(b2),ses(1).time.bins))
hold on
semilogx(ses(1).time.bins,histc(ses(1).time.bins(b1),ses(1).time.bins))

%% Works but looks bad

% plot phase portraits for each axoaxo groups
close all

k = gaussian2Dfilter([20 10],[3 2]);
trans_phase = cellfun(@(a) nanzscore(nanconvn(a,k),[],2),{ph_map(1:end).trans_phase},'uni',0);
trans_phase = cellArrayTo3D(trans_phase);

trans_phase(trans_phase<=0) = 1e-10;
trans_phase = (trans_phase);

k = gaussian2Dfilter([20 10],[1 1]);


ph_rate1 = cellfun(@(a) (nanconvn(a,k)),ph_rate,'uni',0);
ph_rate1 = cellArrayTo3D(ph_rate1);



figure
ax = tight_subplot(3,2);


axes(ax(1))


imagesc(ph_bin,[],(nanmean((ph_rate1(:,:,kp(idx==1))),3)))

hold on
plot(ph_bin,10+cos(ph_bin)*10,'w')
set(gca,'ytick',1:10:length(freq),'yticklabel',round(freq(1:10:end)),'xticklabel',[])
set(gca,'ydir','normal')



axes(ax(2))


imagesc(ph_bin,[],(nanmean((ph_rate1(:,:,kp(idx==2))),3)))

hold on
plot(ph_bin,10+cos(ph_bin)*10,'w')
set(gca,'ytick',1:10:length(freq),'yticklabel',round(freq(1:10:end)),'xticklabel',[])
set(gca,'ydir','normal')


axes(ax(3))


pre_kp = pre_idx(ismember(post_idx(:,[1 4]) ,axax(idx==1,:),'rows'),:);


imagesc(ph_bin,[],(nanmean(ph_rate1(:,:,ismember(cid,pre_kp(:,[1:3]),'rows')),3)))




hold on
plot(ph_bin,10+cos(ph_bin)*10,'w')
set(gca,'ytick',1:10:length(freq),'yticklabel',round(freq(1:10:end)),'xticklabel',[])
set(gca,'ydir','normal')


axes(ax(4))


pre_kp = pre_idx(ismember(post_idx(:,[1 4]) ,axax(idx==2,:),'rows'),:);


imagesc(ph_bin,[],(nanmean(ph_rate1(:,:,ismember(cid,pre_kp(:,[1:3]),'rows')),3)))




hold on
plot(ph_bin,10+cos(ph_bin)*10,'w')
set(gca,'ytick',1:10:length(freq),'yticklabel',round(freq(1:10:end)),'xticklabel',[])
set(gca,'ydir','normal')




axes(ax(5))
kp_con = (ismember(mono_con(:,[1 3]) ,axax(idx==1,:),'rows'));


imagesc(ph_bin,[],(nanmean((trans_phase(:,:,kp_con)),3)))

hold on
plot(ph_bin,10+cos(ph_bin)*10,'w')
set(gca,'ytick',1:10:length(freq),'yticklabel',round(freq(1:10:end)),'xticklabel',[])
set(gca,'ydir','normal')
kp_con = (ismember(mono_con(:,[1 3]) ,axax(idx==2,:),'rows'));

axes(ax(6))
imagesc(ph_bin,[],(nanmean((trans_phase(:,:,kp_con)),3)))

hold on
plot(ph_bin,10+cos(ph_bin)*10,'w')
set(gca,'ytick',1:10:length(freq),'yticklabel',round(freq(1:10:end)),'xticklabel',[])
set(gca,'ydir','normal')


%% Error


%plot mean ripple coupling for presynaptic cells for each axoaxo type
figure

pre_kp = pre_idx(ismember(post_idx(:,[1 4]) ,axax(idx==1,:),'rows'),:);
pre_kp1 = ismember(cid,pre_kp(:,[1:3]),'rows');
plot(-1000:1000,nanmean((rip_ccg(pre_kp1,:))),'k')

hold on

pre_kp = pre_idx(ismember(post_idx(:,[1 4]) ,axax(idx==2,:),'rows'),:);
pre_kp2 = ismember(cid,pre_kp(:,[1:3]),'rows');
plot(-1000:1000,nanmean((rip_ccg(pre_kp2,:))),'r')






%% Error
%plot rate x power corrs

theta_pow = max(amp_cor1(:,35:45),[],2);

kp = find(p<.01 & ratemod>1);
amp_cor1 = cell2mat(amp_cor(2:end)');
close all
figure
ix=1;
ph_bin = linspace(-pi,pi,16);
k = gaussian2Dfilter([10 10],[.5 .5]);


ax = tight_subplot(6,6);
for i = kp'
    axes(ax(ix))
     semilogx(freq(1:end-1),amp_cor1((i),:))

    title(num2str(cid((i),:)))
    ix = ix+1;
    ylim([-.1 .3])
end
%% Works partially then has error

%plot depth for pre and post for each group

bins = -200:30:200;
close all

hold on
plot(bins,histc(celdep(pre_kp1),bins)/sum(pre_kp1),'k')
plot(bins,histc(celdep(pre_kp2),bins)/sum(pre_kp2),'r')

bins = -200:40:200;

figure
hold on
 plot(bins,histc(celdep(kp(idx==1)),bins)/sum(idx==1),'k')
 plot(bins,histc(celdep(kp(idx==2)),bins)/sum(idx==2),'r')



%% Works partially then has error

 %plot depth x phase locking/ripple mode/theta corr
 close all
 figure
 plot(ripmod(kp(idx==1)),celdep(kp(idx==1)),'o','color','k')
 hold on
 plot(ripmod(kp(idx==2)),celdep(kp(idx==2)),'x','color','r')
 xlabel('Ripple mod')
 ylabel('Depth')
 figure
 plot(theta_pow(kp(idx==1)),celdep(kp(idx==1)),'o','color','k')
 hold on
 plot(theta_pow(kp(idx==2)),celdep(kp(idx==2)),'x','color','r')

  xlabel('Theta corr')
 ylabel('Depth')
 figure
 plot(ph_pref_theta(kp(idx==1)),celdep(kp(idx==1)),'o','color','k')
 hold on
 plot(ph_pref_theta(kp(idx==2)),celdep(kp(idx==2)),'x','color','r')
 xlabel('Theta pref')
 ylabel('Depth')

 %% Works

 %plot connectivity maps
  pre_kp = pre_idx(ismember(post_idx(:,[1 4]) ,axax(idx==1,:),'rows'),:);

  figure
   hold on
  for i = 1:length(axax)

      pre_kp = pre_idx(ismember(post_idx(:,[1 4]) ,axax(i,:),'rows'),:);
       rr(1) = rand(1)/2+1;
          if idx(i) ==1
        plot(rr(1),celdep(kp(i)),'o','color','k')
          else
         plot(rr(1),celdep(kp(i)),'o','color','r')
          end
      for j = 1:size(pre_kp,1)
          rr(2) = rand(1)/2;
       pre_kp1 = ismember(cid,pre_kp(j,[1:3]),'rows');

          if idx(i) ==1
 plot([0 0]+rr,[celdep(kp(i)),celdep(pre_kp1)],'k')



 plot(rr(2),celdep(pre_kp1),'^','color','k')
          else
               plot([0 0]+rr,[celdep(kp(i)),celdep(pre_kp1)],'r')



 plot(rr(2),celdep(pre_kp1),'^','color','r')
          end

      end
  end



