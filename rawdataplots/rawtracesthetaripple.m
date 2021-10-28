%% Unc5b FR x Other


% load([basename '.cell_metrics.cellinfo.mat']) % file cellexplorer

% cellTypeLabels = cell_metrics.putativeCellType;
% for iL = 1:length(cellTypeLabels)
%    
%     NarrowINT(iL) = strcmpi(cellTypeLabels{iL},'Narrow Interneuron');
%     WideINT(iL) = strcmpi(cellTypeLabels{iL},'Wide Interneuron');
%     Pyr(iL) = strcmpi(cellTypeLabels{iL},'Pyramidal Cell');
% end
% 
% nINTInd = find(NarrowINT);
% AACInd = [3 13 61];
% wINTInd = find(WideINT);
% PyrInd = find(Pyr);

[pyrs, ints, aacs] = splitCellTypes(basepath);
%% Example theta trace + rasters AACs underneath
% if strcmpi(basename,'u19_200313_155505')
  basename =  'mouse1_180501a'
opts.rippleChan     = 27% 0 was specified, changed it after looking at the dat file.%57;
% end
chanRip = 27;
[ripples] = bz_FindRipples(cd,chanRip,'durations',[50 150],...
'thresholds',[1 3], 'EMGThresh', 0.95,'saveMat',true);

lfp = bz_GetLFP(opts.rippleChan);
spikes = bz_LoadPhy

%%
selXLim = [3010 3020];
figure,
s1 = subplot(3,1,1);
plot(lfp.timestamps,lfp.data)
xlim(selXLim)

% get and plot theta
s2 = subplot(3,1,2);
lfp_theta = BandpassFilter(double(lfp.data), 1250 ,[6 8]);
plot(lfp.timestamps,lfp_theta*0.195)
xlim([selXLim])
xlabel('time (s)')
ylabel('amplitude (uV)')

% now plote AAC spike times in it
s3 =subplot(3,1,3);
% startUnder = -800%min(lfp_theta);
% plotShift = 0;

for iSel = 1:length(aacs)
    selAAC = aacs(iSel);
    x = spikes.times{selAAC};
    totSpkTims = length(spikes.times{selAAC});
    y = [repmat(iSel-1,totSpkTims,1),repmat(iSel,totSpkTims,1)];
    hold on
     
    tx = [x.';x.';nan(1,length(x))];
    ty = [y(:,1).';y(:,2).';nan(1,length(x))];
    plot(tx(:),ty(:))
    xlim([selXLim])
xlabel('time (s)')
ylabel('AAC')
end

for iThEpoch = 1:100:lfp.timestamps(end)
    widthTE = 1;
    thEpoch = [iThEpoch-widthTE iThEpoch+widthTE];
    s1.XLim=thEpoch;
    s2.XLim=thEpoch;
    s3.XLim=thEpoch;
    saveDir = ['E:\Dropbox\PD_Hpc\Progress\AAC\examples_theta_rip_' basename '\thetanew\'];
    sD = dir(saveDir);
    if isempty(sD)
        mkdir(saveDir)
    end
          print(gcf,[saveDir 'theta_example_' num2str(iThEpoch) '.pdf'], '-dpdf','-bestfit')

   end


  %%
   load([basename '.ripples.events.mat'])
  for iRip = 1:length(ripples.peaks)
      
      selRipple = iRip;
      timeAroundripple = .2;
      ripXLim = [ripples.peaks(selRipple)-timeAroundripple ripples.peaks(selRipple)+timeAroundripple];
      
      figure,
      subplot(3,1,1)
      plot(lfp.timestamps,lfp.data)
      xlim(ripXLim)
      
      % get and plot theta
      subplot(3,1,2)
      % get and plot ripple
      lfp_ripple = BandpassFilter(double(lfp.data), 1250 ,[100 250]);
      plot(lfp.timestamps,lfp_ripple*0.195)
      
      xlim(ripXLim)
      xlabel('time s')
ylabel('amplitude uV')
      subplot(3,1,3)
      
      for iSel = 1:length(aacs)
          selAAC = aacs(iSel);
          x = spikes.times{selAAC};
          totSpkTims = length(spikes.times{selAAC});
          y = [repmat(iSel-1,totSpkTims,1),repmat(iSel,totSpkTims,1)];
          hold on
          
          tx = [x.';x.';nan(1,length(x))];
          ty = [y(:,1).';y(:,2).';nan(1,length(x))];
          plot(tx(:),ty(:))
          xlim(ripXLim)
          xlabel('time s')
ylabel('AAC')
      end
      saveDir = ['E:\Dropbox\PD_Hpc\Progress\AAC\examples_theta_rip_' basename '\ripnew\'];
         sD = dir(saveDir);
    if isempty(sD)
        mkdir(saveDir)
    end
      print(gcf,[saveDir 'ripple_example_' num2str(selRipple) '.pdf'], '-dpdf','-bestfit')
      close
  end