%% Workflow_axax_help_for_lianne_u19_200313_155505(2)
%

%% get rip idxs
rip_idx = Find_closest_idx(ripples.peaks, lfp.timestamps)

%% load lfp
chans = [58 37 19 4]; % best rip chans from lianne
lfp = bz_GetLFP('all');
%% average rip traces

win = 250;
rip_pul_avg = Running_sum_avg(lfp.data, win, rip_idx(aRP));
rip_nopul_avg = Running_sum_avg(lfp.data, win, rip_idx(~aRP));
%% Voltage conversion
basename = bz_BasenameFromBasepath(cd);
num_channels = length(lfp.channels);
fileinfo        = dir([basename '_analogin.dat']);
num_samples     = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes


fid = fopen([basename '_analogin.dat'], 'r');
v   = fread(fid, [num_channels, num_samples], 'uint16');
fclose(fid);
v   = v * 0.000050354; % convert to volts, intan conversion factor 
%% plot rip averages

figure

for i = 1:length(lfp.channels)
    subplot(8,8,i)
   
    plot(rip_pul_avg{i}*v,'k')
    hold on
    plot(rip_nopul_avg{2}*v,'r')
    title([num2str(lfp.channels(i))])
    ylabel(['mV'])
    xlabel(['Time (s)'])
    set(gca,'XTick',[0 win/2 win+1 win*1.5 win*2+1])
    set(gca,'XTickLabel',[-win/1250 -win/1250/2 0 win/1250/2 win/1250])
    set(gca,'XLim',[0 win*2+1])
    set(gca,'YLim', [-0.25 0.1])
%     set(gca,'YLim',[-.14 .07])
end

%% Getting lianne's cells to match sam's

% Works!!
mono_res = bz_GetMonoSynapticallyConnected(cd, 'saveMat', true, 'plot', false);

%%
ShortTermPlasticity_axax_Kaiser_edit(cd)




































