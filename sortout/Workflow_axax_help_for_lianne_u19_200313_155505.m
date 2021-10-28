%% Workflow_axax_help_for_lianne_u19_200313_155505

lfp = bz_GetLFP('all') %load in LFP file

x = ismember(pulseEpochs, double(lfp.timestamps)) %verify all times match from the recording
sum(x) % make sure numbers match between x and pulseEpochs

%% need to find the closest value that matches

str_idx = []; % reserve/ empty

for i = 1:length(pulseEpochs) % Loop through pulses
    [~,min_idx] = min(abs(pulseEpochs(i,1) - lfp.timestamps)); % Find closest matching time
    str_idx(i,1) = min_idx; % store
end

%% filter lfp for chans
chans = [58 37 19 4]; % best rip chans from lianne
BP_rng = [120, 200];
LP_rng = [250];

best_chans = lfp.data(:,chans);
for i = 1:size(best_chans,2)
    BP_LFP(:,i) = bandpass(double(best_chans(:,i)),[120, 200],1250); % bandpass filter 
end

% lowpass filter for SWs
for i = 1:size(best_chans,2)
    LP_LFP(:,i) = lowpass(double(best_chans(:,i)),LP_rng,1250); % lowpass filter 
end


%% pull data for average window

win = 625; % in samples equaling 500 ms
pul_data_rip_BP = []; % reserve/ empty

for ii = 1:length(chans)
    for i = 1:length(str_idx); % Loop through indices
        pul_data_rip_BP{ii}(i,:) = double(BP_LFP(str_idx(i) - win:str_idx(i) + win,ii)); % Collect data from the window
    end
end

for ii = 1:length(chans)
    for i = 1:length(str_idx); % Loop through indices
        pul_data_sw_BP{ii}(i,:) = double(LP_LFP(str_idx(i) - win:str_idx(i) + win,ii)); % Collect data from the window
    end
end

%% Intan voltage conversion

basename = bz_BasenameFromBasepath(cd);
num_channels = length(lfp.channels);
fileinfo        = dir([basename '_analogin.dat']);
num_samples     = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes


fid = fopen([basename '_analogin.dat'], 'r');
v   = fread(fid, [num_channels, num_samples], 'uint16');
fclose(fid);
v   = v * 0.000050354; % convert to volts, intan conversion factor 

%% Plot Rip_BP



figure

for i = 1:length(pul_data_rip_BP)
    subplot(2,2,i)
    for ii = 1:size(pul_data_rip_BP{i},1)
        plot(pul_data_rip_BP{i}(ii,:),'color', [0,0.7,0.9,.5])
        hold on
    end
    line([win+1 win+1], [-4000 4000], 'color', 'green')
    line([win+1+(win*.3*2) win+1+(win*.3*2)], [-4000 4000], 'color', 'red')
    plot(mean(pul_data_rip_BP{i},1),'k')
    title([num2str(chans(i))])
    ylabel(['NA'])
    xlabel(['Time (s)'])
    set(gca,'XTick',[0 win/2 win+1 win*1.5 win*2+1])
    set(gca,'XTickLabel',[-win/1250 -win/1250/2 0 win/1250/2 win/1250])
    set(gca,'XLim',[0 win*2+1])
%     set(gca,'YLim',[-.14 .07])
end
%% SW_BP
figure

for i = 1:length(pul_data_sw_BP)
    subplot(2,2,i)
    for ii = 1:size(pul_data_sw_BP{i},1)
        plot(pul_data_sw_BP{i}(ii,:),'color', [0,0.7,0.9,.5])
        hold on
    end
    line([win+1 win+1], [-4000 4000], 'color', 'green')
    line([win+1+(win*.3*2) win+1+(win*.3*2)], [-4000 4000], 'color', 'red')
    plot(mean(pul_data_sw_BP{i},1),'k')
    title([num2str(chans(i))])
    ylabel(['NA'])
    xlabel(['Time (s)'])
    set(gca,'XTick',[0 win/2 win+1 win*1.5 win*2+1])
    set(gca,'XTickLabel',[-win/1250 -win/1250/2 0 win/1250/2 win/1250])
    set(gca,'XLim',[0 win*2+1])
%     set(gca,'YLim',[-.14 .07])
end


%% Load in specific channels for lfp
chans = [57 36 18 3];

lfp = bz_GetLFP(chans); % When loading in these channels use base-0 indexing

%% Load all channels
lfp = bz_GetLFP([0:15]);
%% Wave spec

temp_lfp = lfp;
temp_rip_ws = [];



for i = 1:length(str_idx) % loop through idices

    for ii = 1:size(lfp.data,2) % loop through data
        tic
        temp_lfp.data = double(lfp.data(str_idx(i)-250:str_idx(i)+550,ii)./lfp.data(str_idx(i)-250-801:str_idx(i)-250-1,ii)); % collect data of specific window, normalize to window pre-stim and store in temp lfp
        temp_lfp.timestamps = lfp.timestamps(str_idx(i)-250:str_idx(i)+550); % collect corresponding times from the window
        ws{i,ii} = bz_WaveSpec(temp_lfp, 'frange', [1 250], 'nfreqs', [250],... % run wavespec and store output
            'space', 'lin', 'samplingRate', [1250]);
        toc
    end
    


end

%% average all the output from the wavespec

for i = 1:size(ws,2) % loop through channels
    run_sum = zeros;
    for ii = 1:length(ws)
        run_sum = run_sum+abs(ws{ii,i}.data);
    end
    mean_ws{1,i} = run_sum/ii;
end


%% Plot WS
str = 250;
stp = 550;

figure
for i = 1:length(mean_ws)
    subplot(4,4,i)
    imagesc(mean_ws{i}')
    title([lfp.channels(i)]) % '' num2str(volt(ivolt))])
    xlabel('Time(ms)')
    ylabel('Frequency(Hz)')
    % caxis([0 5])
    set(gca, 'YDir', 'normal')
    set(gca, 'YTick', [1:50:250])
    set(gca, 'YTickLabel', {1:50:250})
    set(gca, 'XTick', [1 201 401 601 801])
    set(gca, 'XTickLabel', {-400 -200 0 200 400})
    colorbar
    caxis([0 300])
    hold on
    line([str+1 str+1], [1 25], 'color', 'green')
    line([stp+1 stp+1], [1 25], 'color', 'red')
end





