

%%

%make pulses


%make an event file  with basename.evt.sti

dirName  = pwd;
sl = regexp(dirName,'/');
basename = dirName(sl(end)+1:end);
analogfils_dat = [basename '_analogin.dat'];
%define experiment start and stop based off of pulses -  pulse_on pulse_off

%define the basename.evt.sti with pulse_on/pulse_off
system(['neuroscope ' analogfils_dat])



%%
%%get pulse times

dirName  = pwd;
sl = regexp(dirName,'/');
basename = dirName(sl(end)+1:end);
analogfils_dat = [basename '_analogin.dat'];


filename = getAllExtFiles(pwd,'sti',0);
events = LoadEvents(filename{1});
MakeEventTime([dirName '/' analogfils_dat],events,'ch',[2],'thres1',500,'thres2',12500); %(1:# is channels 1-#, or [2 4] if only 2 and 4)



%%


%make LFP
bz_LFPfromDat(pwd);

%%

% get ripples, choose a good channel (base 0)
[ripples] = bz_FindRipples(pwd,0,'thresholds',[2 4],'basename','mouse2_180523_lfp','frequency',1250,'nChannels',2)



%%
%load spikes


spikes = bz_GetSpikes('getWaveforms',false,'noprompts',true);
basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);

ripFil = [basepath '/' basename '.evt.rip'];
LFPFil = [basepath '/' basename '.lfp'];
rip_evs = LoadEvents(ripFil);

rip_ts = rip_evs.time(cellfun(@any,regexp(rip_evs.description,'peak')));

rip_st = rip_evs.time(cellfun(@any,regexp(rip_evs.description,'start')));

rip_end = rip_evs.time(cellfun(@any,regexp(rip_evs.description,'stop')));

stims = LoadEvents([basename '.evt.ait']);
on = stims.time(cellfun(@any,regexp(stims.description,'pulse_on')));
off = stims.time(cellfun(@any,regexp(stims.description,'pulse_off')));

basename = bz_BasenameFromBasepath(pwd);



%%
% load clu fil, if you have multiple shanks you will need to load separate
% files per shank
 basename = bz_BasenameFromBasepath(pwd);
[xml, rxml] = LoadXml([basename  '.xml']);
nShank = length(xml.SpkGrps);
clu = [];
for i  = 1:nShank
    
clu{i} = LoadClu([basename '.clu.' num2str(i)]);
end
%%
% make 
loner = rip_st(diff([0;rip_st])>.5);
rip_peth = [];
for i = 1:length(spikes.times)
rip_peth(i,:) = CrossCorr(loner,spikes.times{i},.01,100);
end

 rip_peth = (rip_peth/length(loner)*100);
%%

 %2, 29
stim_peth = [];
for i = 1:length(spikes.times)
stim_peth(i,:) = CrossCorr(on,spikes.times{i},.01,100);
end


 stim_peth = (stim_peth/length(on)*100);
%%
 
 acg =[];rate =[];
 for i = 1:length(spikes.times)
 [status, interval] = InIntervals(spikes.times{i},[on off]);

 spks = spikes.times{i}(~status);
 
 ok = CrossCorr(spks,spks,.001,200);
 
 acg(i,:) = ok/sum(~status);
 acg(i,101) = nan;
 
 rate(i) = sum(~status)/spikes.spindices(end,1);
 end
 




 %%
 
 
fname = [basename '.spk.1'];
waveforms = LoadSpikeWaveforms(fname,32,48);
waveforms = permute(waveforms,[2 3 1]);
%%
%%
mono_res = MonoSynClick([],'plot',false);
%%
 
%get mean waveform

pos = [];ix = 1;

uclu = setdiff(unique(clu),0)';
for i = uclu
    
wv = nanmean(waveforms(:,:,clu==i),3);
[~,b] = min(min(wv,[],2));

pos(ix,:) = [xx(b) yy(b)];
ix = ix+1;
end
 
 %%
 

 close all
 int = [1 2 3 4 5 7 8 9 19];
 ax = [1 2 3 7 8 9];
 pyr = setdiff(1:length(uclu),[int]);
 
 kp = ismember(mono_res.sig_con(:,2),int) & ismember(mono_res.sig_con(:,1),pyr);
 
 sig_con = mono_res.sig_con(kp,:);
 
 figure
 pos1 =[];
 for i = 1:length(pos)
 pos1(i,:) = pos(i,:)+rand(1)/2;
 end
 for i = 1:length(sig_con)
     
     if ismember(sig_con(i,2),ax)
          plot([pos1(sig_con(i,1),1) pos1(sig_con(i,2),1)]',[pos1(sig_con(i,1),2) pos1(sig_con(i,2),2)]','color',[.5 .5 1])
     else
         
     plot([pos1(sig_con(i,1),1) pos1(sig_con(i,2),1)]',[pos1(sig_con(i,1),2) pos1(sig_con(i,2),2)]','color',[1 .5 .5])
     end
     hold on
     
 end
 
 for i = 1:length(pyr)
     plot(pos1(pyr(i),1),pos1(pyr(i),2),'^','color','k')
 end
 
 
 
 for i = 1:length(int)
     plot(pos1(int(i),1),pos1(int(i),2),'o','color','r')
 end
 
 
 
 for i = 1:length(ax)
     plot(pos1(ax(i),1),pos1(ax(i),2),'x','color','b')
 end
 
 
 axis square
 axis off
 set(gca,'ydir','reverse')
 
 %%
 
 

%%
% 1,2,3 6,7,8
cel = 4;

xx = [2 repmat(1:3,1,10) 2];
yy = [1 upSample(2:11,3)' 12];
figure 
 ax = tight_subplot(1,4);
 
 
 axes(ax(1))
 
 
 for i = 1:32
 plot(100*xx(i)+(1:48), squeeze(nanmean(waveforms(i,:,clu==cel+1),3))-yy(i)*700,'k')
 hold on
 end
 axis off
 
 
 axes(ax(2))
 bar(-100:100,acg(cel,:),'k')
xlim([-100 100])



 
 axes(ax(3))
 bar(-50:50,stim_peth(cel,:),'k')
xlim([-50 50])
ylim([0 150])

 axes(ax(4))
 bar(-50:50,rip_peth(cel,:),'k')
xlim([-50 50])
ylim([0 20])


mtit(num2str([spikes.shankID(cel) spikes.cluID(cel)]))


%%


%  %%
%   [status, interval] = InIntervals(spikes.times{16},[on off]);
% 
%  spks = spikes.times{16}(~status);
%  
%  ok = CrossCorr(spks,spks,.001,200);
%  ok = ok([1:100 102:end]);
%  
%  
%  mono_res.sig_con;
%  
%  %%
%  
%  
% preax =CrossCorr(spikes.times{35},spikes.times{20},.001,100);
% prenonax =CrossCorr(spikes.times{35},spikes.times{6},.001,100);


%%
% ph = [];
% ts = [];
% 
% for i = 1:length(rip_ts)
%     
%     
%     st = rip_ts(i);
%     pr = round(1250*(st - rip_st(i)));
%      po = round(1250*(rip_end(i) - st));
%     data = LoadBinary(LFPFil,'frequency',1250,'channels',14,'start',st-1,'duration',2,'nchannels',32);
%     
%     dataf = BandpassFilter(double(data),1250,[80 220]);
%     pht = InstPhase(dataf);
% % [wavespec] = bz_WaveSpec(double(data),'samplingRate',1250,'frange',[100 120], 'nfreqs' ,1);
% %get phase at max power
%  
% %InstPhase = atan2(imag(wavespec.data), real(wavespec.data))';
% %[a,b] = max(abs(wavespec.data)');
% 
% %idx = sub2ind(size(InstPhase),b,1:size(InstPhase,2));
% %InstPhase = InstPhase(1,:);
%  
% ph =  [ph pht((-pr:po)+1251)'];
% ts = [ts st+((-pr:po)+1251)/1250];
% end
% [a,b] = histc(spikes.times{20},ts);
% b = b(a(b>0)==1);
% 
% bin_spk = histc(ph(b),-pi:pi/10:pi);
% bin_ph = histc(ph,-pi:pi/10:pi);
% rip_rate = (bin_spk./bin_ph)*1250;
% 
% bar(-pi:pi/10:pi,rip_rate)

 
 