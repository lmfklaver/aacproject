% a = [0 1.2 1.3 1.4 1.5 2 3 5.11 5.2 5.3 16 17.1 17.2 19]

% Exclude pulses
pulseEpochs = optoStim.timestamps;

% Get spikes in or out pulse epoch
[status_pulse ,~ , ~ ] = cellfun(@(a) InIntervals(a,pulseEpochs), spikes.times,'UniformOutput', false);

for iUnit = 1:length(spikes.times)
    %         spkTimIN{iUnit}   = spikes.times{iUnit}(status_pulse{iUnit});
    spkTimNoPulse{iUnit}   = spikes.times{iUnit}(~status_pulse{iUnit});
end
%%
iUnit = 1;
a = spkTimNoPulse{iUnit};
%change to gd_eps
bursty = [];
    for jj = 2 : length(a) - 1
        bursty(jj) =  any(diff(a(jj-1 : jj + 1)) < 0.006); %less than 6ms
    end
 myCalc = length(find(bursty > 0))/length(bursty) 
 
%  spksInBurst = a([logical(bursty)])
 CE =cell_metrics.burstIndex_Mizuseki2012(iUnit)