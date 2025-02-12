
launchDirNforAACSessions

sessions = [24 25];

% % % % % % % % % % % % % % % % % % % % % % % %
% % % Get cumulative values for SW and MUA
% % %
% % % % % % % % % % % % % % % % % % % % % % % %

cumul_rateSpkperRipON=[];
cumul_rateSpkperRipOFF=[];
cumul_sharpwavepeakUvON=[];
cumul_sharpwavepeaknormON=[];
cumul_sharpwavepeakZON=[];
cumul_sharpwavepeakPeakZON=[];
cumul_sharpwavepeakUvOFF=[];
cumul_sharpwavepeaknormOFF=[];
cumul_sharpwavepeakZOFF=[];
cumul_sharpwavepeakPeakZOFF=[];    
cumul_RipPOWON=[];
cumul_RipPOWOFF=[];
    
    
for iSess = sessions
    cd(dirN{iSess})
    basepath    = cd;
    basename    = bz_BasenameFromBasepath(basepath);
    load([basename '_celltypes.mat']);
    load([basename '.ripspikes.allripinstim.analysis.mat'])
    load([basename '.spikes.cellinfo.mat'])
    load([basename '.ripples.events.mat'])
    pyrspikes=[];
    rateSpkperRipON=[];
    rateSpkperRipOFF=[];
    for i=1:length(pyrs)
    pyrspikes=[pyrspikes;spikes.times{i}];
    end
    pyrspikes=sortrows(pyrspikes);
    
    [statusON,intervalON,indexON]=InIntervals(pyrspikes,ripspikes.ONrips.timestamps);
    [statusOFF,intervalOFF,indexOFF]=InIntervals(pyrspikes,ripspikes.OFFrips.timestamps);
    rateSpkperRipON=zeros(1,size(ripspikes.ONrips.timestamps,1));
    for iInterval = unique(intervalON(intervalON~=0))';
        numSpkON= sum(length(find((intervalON==iInterval))));
        rateSpkperRipON(iInterval)=(numSpkON/(ripspikes.ONrips.timestamps(iInterval,2)-ripspikes.ONrips.timestamps(iInterval,1)))/size(ripspikes.spikesRipNum,1);
    end
    rateSpkperRipOFF=zeros(1,size(ripspikes.OFFrips.timestamps,1));
    for iInterval = unique(intervalOFF(intervalOFF~=0))';
        numSpkOFF= sum(length(find((intervalOFF==iInterval))));
        rateSpkperRipOFF(iInterval)=(numSpkOFF/(ripspikes.OFFrips.timestamps(iInterval,2)-ripspikes.OFFrips.timestamps(iInterval,1)))/size(ripspikes.spikesRipNum,1);
    end
    
    cumul_rateSpkperRipON=[cumul_rateSpkperRipON rateSpkperRipON];
    cumul_rateSpkperRipOFF=[cumul_rateSpkperRipOFF rateSpkperRipOFF];
    cumul_RipPOWON=[cumul_RipPOWON ripspikes.ONrips.peakNormedPower'];
    cumul_RipPOWOFF=[cumul_RipPOWOFF ripspikes.OFFrips.peakNormedPower'];
    cumul_sharpwavepeakUvON=[cumul_sharpwavepeakUvON ripspikes.ONrips.sharpwavepeakUv'];
    cumul_sharpwavepeaknormON=[cumul_sharpwavepeaknormON ripspikes.ONrips.sharpwavepeaknorm'];
    cumul_sharpwavepeakZON=[cumul_sharpwavepeakZON ripspikes.ONrips.sharpwavepeakZ'];
    cumul_sharpwavepeakPeakZON=[cumul_sharpwavepeakPeakZON ripspikes.ONrips.SWpeakZScore'];
    cumul_sharpwavepeakUvOFF=[cumul_sharpwavepeakUvOFF ripspikes.OFFrips.sharpwavepeakUv'];
    cumul_sharpwavepeaknormOFF=[cumul_sharpwavepeaknormOFF ripspikes.OFFrips.sharpwavepeaknorm'];
    cumul_sharpwavepeakZOFF=[cumul_sharpwavepeakZOFF ripspikes.OFFrips.sharpwavepeakZ'];
    cumul_sharpwavepeakPeakZOFF=[cumul_sharpwavepeakPeakZOFF ripspikes.OFFrips.SWpeakZScore'];

    end

