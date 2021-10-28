% launchDirNforAACSessions

for iSess = 1:2%5;
    cd(dirN{iSess})
    basepath = cd;
basename = bz_BasenameFromBasepath(cd);
% load([basename '.cell_metrics.cellinfo.mat'])
% load([basename '.spikes.cellinfo.mat'])

%% Ripple masking Workflow
load([basename '.ripples.events.mat'])
load([basename '.spikes.cellinfo.mat'])

%%
% I want to collect ripples that have a peak that is more than 5 seconds
% from the end of the previous ripple and more than 5 seconds to the next
% ripple

 rip_times = [];
% place rip times and peak all together
rip_times(:,1) = ripples.peaks;
rip_times(:,2:3) = ripples.timestamps;

lfp = bz_GetLFP(ripples.detectorinfo.detectionparms.channel);

rip_ind = find(ismember(lfp.timestamps,rip_times(:,1)));
%%
rip_plus = [];
rip_minus = [];
for i = 1:length(rip_times)
    if i == 1 || i == length(rip_times)
        rip_plus = [rip_plus; nan];
        rip_minus = [rip_minus; nan];
    else
        rip_plus = [rip_plus; abs(rip_times(i,1) - rip_times(i+1,2))];
        rip_minus = [rip_minus; abs(rip_times(i,1) - rip_times(i-1,3))];%abs(rip_times(i,1) - rip_times(i+1,2)) > 5 &&  abs(rip_times(i,1) - rip_times(i-1,3)) > 5
%         ripind = [ripind;i];
    end
end

rip_plus_samps = abs(round(rip_plus * lfp.samplingRate));
rip_minus_samps = abs(round(rip_minus * lfp.samplingRate));
for i = 1:length(rip_plus_samps)
   
    if rip_minus_samps(i) > 6250
        rip_minus_samps(i) = 6250;
       
    end
    if rip_plus_samps(i) > 6250
        rip_plus_samps(i) = 6250;
    end
end

%% firing rate for each bin of time for each cell during each ripple

% 10ms bins need to be variable since the windows for ripples are also
% variable


%for ii = 1:length(spikes.times)

%%
Cell_Spk_bins = [];
%%%%%
[pyrs, ints, aacs] = splitCellTypes(basepath);
%%%%%
AACCount = 0;
for i = aacs%1:spikes.numcells
    AACCount =AACCount+1;
    tic
    Spk_bins = nan(length(rip_ind), 1000 +1);
    for ii = 1:length(rip_ind)
        if isnan(rip_minus(ii)) || isnan(rip_plus(ii))
            continue
        end
        temp_bins = 0:.01:rip_minus(ii);
        if length(temp_bins)>500
            temp_bins = temp_bins(1:500);
        end
        temp_bins = flip(temp_bins,2);
        temp_pre_bins = rip_times(ii,1)-temp_bins;
       
        temp_bins = .01:.01:rip_plus(ii);
        if length(temp_bins)>500
            temp_bins = temp_bins(1:500);
        end
       
        temp_post_bins = rip_times(ii,1)+temp_bins;
       
        temp_bins = [temp_pre_bins temp_post_bins];
       
        temp_ints = [];
       
        temp_ints(:,1) = temp_bins(1:end-1)';
        temp_ints(:,2) = temp_bins(2:end)';

        temp = nan(1,length(temp_ints));
        [status,interval,index] = InIntervals(spikes.times{i},temp_ints);
        temp(nonzeros(interval)) = nonzeros(index);
       
        Spk_bins(ii,(((size(Spk_bins,2)-1)/2 +1)-length(temp_pre_bins))+1:(((size(Spk_bins,2)-1)/2 +1)+length(temp_post_bins))-1) = temp;
    end
    Cell_Spk_bins(:,:,AACCount) = Spk_bins;
    toc
end
 M{iSess} = Cell_Spk_bins;
end

% % %%
% % figure
% % hold on
% % for i = 1:length(Cell_Spk_bins)
% %     plot(nanmean(Cell_Spk_bins{i},1))
% % end
