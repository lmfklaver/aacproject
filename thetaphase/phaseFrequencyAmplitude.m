%% phaseFrequencyAmplitude
%
%  MAYBE I CAN JUST DO THIS WITH THE INPUTS IN GETPHASEMAP?!?! CHECK IT OUT
%
%
%
%
% determine frequency theta ([4 10])
% get phases of all spikes to each of those frequencies

freq = [4:1:10];

% for i = 1%:nFreq-1

    %filter LFP
    filter.data = BandpassFilter(double(d), 1250 ,[freq(1) freq(end)]); % Filter the data for the specified frequency made above
    
    %get inst phase and amplitude
    filtered.phase  = InstPhase(filter.data ); % inst phase using the hilbert transform
    filtered.amp    = InstAmplitude(filter.data ); % inst amp using hilbert transorm
    filtered.timestamps     = ts ; 
    filtered.sampleRate   = xml.lfpSampleRate;
    filtered.filterparms.passband = [freq(1) freq(end+1)]; % record passband range of the filter for that loop
    spikes1.times{1} = sort(cell2mat(spikes.times')); 
    % puts all cells from the recording into one variable to calculate MUA
    % for determining powerbins to keep (Sam was not sure if this is a good
    % idea)
    
[PowerPhaseRatemap,spikebinIDs] = bz_PowerPhaseRatemap_sam(spikes,filtered,filtered,'powernorm','none');

% end
