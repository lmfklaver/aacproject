for iSess = 1%:length(dirN)
    
    cd(dirN{iSess})
    
    basepath =cd;
    basename = bz_BasenameFromBasepath(basepath);
    
    
    load([basename '.STP.mat'])
    load([basename '.optoStim.manipulation.mat'])
    load([basename '.run.states.mat'])
    
    
    %% get correct LFP channel (one from ripple)
    
    % get ripple
    fils  = getAllExtFiles(basepath,'rip',1);
    rip = LoadEvents(fils{1});
    
    % pulls the channel from the ripples and loads the xml file
    ch = str2num(rip.description{1}(regexp(rip.description{1},'[0-9]')));
    xmln = [basepath '/' basename '.xml'];
    fname = [basepath '/' basename '.lfp'];
    xml = LoadXml(xmln);
    
    % load lfp
    d =[];ts =[];
    dt = LoadBinary(fname,'nchannels',xml.nChannels,'channels',ch+1,'frequency',xml.lfpSampleRate);
    d  = [d;dt];
    ts = [ts ; [0:length(dt)-1]'/1250];
    
      
    %load spikes
    [spikes] = bz_LoadPhy;
    %get spike times OUT pulse
    [status, interval] = cellfun(@(a) InIntervals(a,optoStim.timestamps),spikes.times, 'uni',false);
    
    for iUnit = 1:length(spikes.times)
        spikes.nopulse{iUnit} = spikes.times{iUnit}(~status{iUnit});
    end
    
    [status, interval] = cellfun(@(a) InIntervals(a,run.epochs),spikes.nopulse, 'uni',false);
    
    for iUnit = 1:length(spikes.nopulse)
        spikes.run{iUnit} = spikes.nopulse{iUnit}(status{iUnit});
    end
    
    for iUnit = 1:length(spikes.run)
        
        spkInstPhaseAmp.times = {spikes.run{iUnit}};
        
        %define the number of frequencies with which to look at phase
        %preference
        nfreq   = 25; % was echt 200 of zo LK % set number of frequencies
        freq    = logspace(log10(1),log10(250),nfreq); % set frequencies along a log scale
        ph_bin  = linspace(-pi,pi,32); % set phases into 32 bins
        ph_rate = nan(nfreq-1,32,length(spkInstPhaseAmp.times)); % reserve ph_rate space for the number of frequencies and how they are binned
        %%
        %hardcoded for theta
        freq = 5:8;
        
        %filter LFP
        filtered.data = BandpassFilter(double(d), 1250 ,[freq(1) freq(end)]); % Filter the data for the specified frequency made above
        %get inst phaes and amplitude
        filtered.phase  = InstPhase(filtered.data); % pull the instantaneous phase using the hilbert transform
        filtered.amp    = InstAmplitude(filtered.data); % pull the instantaneous amplitude using hilbert transorm
        filtered.timestamps     = ts ; % set timestamps of the filtered signal
        filtered.samplingRate   = 1250; % sampling rate
        filtered.filterparms.passband   = [freq(1) freq(end)]; % record passband range of the filter for that loop
        
        %these are now spikes in run 5cm/s and outside the pulse
        spkInstPhaseAmp.times{iUnit}     = sort(cell2mat(spkInstPhaseAmp.times')); % overbodig, want maar 1 cell - puts all cells from the recording into one variable I am not sure if this is a good idea
        spkInstPhaseAmp.UID  = spikes.UID(iUnit);
        spkInstPhaseAmp.cluID = spikes.cluID(iUnit);
        
        
        %%
        %Get Power/Phase at each spike
        spkInstPhaseAmp.amp{iUnit} = cellfun(@(X) interp1(filtered.timestamps,filtered.amp,X,'nearest'),...
            spkInstPhaseAmp.times,'uniformoutput',false);
        spkInstPhaseAmp.phase{iUnit} = cellfun(@(X) interp1(filtered.timestamps,filtered.phase,X,'nearest'),...
            spkInstPhaseAmp.times,'uniformoutput',false);
        spkInstPhaseAmp.ph_bin{iUnit} = ph_bin;
%         histoEdges = ph_bin;
%         binnedHisto = histcounts(spkInstPhaseAmp.phase{1},histoEdges);
%         figure 
%         histogram(spkInstPhaseAmp.phase{iUnit}{1})
%         title([basename ' phase - minRunSpeed ' num2str(run.detectorinfo.detectionparms)])
%         savefig([basename '_phase'])
%         
%         figure 
%         histogram(spkInstPhaseAmp.amp{iUnit}{1})
%         title([basename ' amp - minRunSpeed ' num2str(run.detectorinfo.detectionparms)])
%         savefig([basename '_amp'])
        
         save([basename '.phamp.analysis.mat'],'spkInstPhaseAmp')
        % histogram(binnedHisto,histoEdges)
    end
end

    
