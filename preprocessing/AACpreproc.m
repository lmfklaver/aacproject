% unc5b env inhibition
basepath = 'D:\Data\Axoaxonic_Data_Lianne';

sessions = {'u19_200310_135409';...
'u19_200313_120452';...      
'u19_200313_155505';...          
'u21_200305_153604';...
'u21_200309_142534';...          
'u26_200306_172032';...          
% 'u26_200308_144613';...  %trash, needs more preprocessing
};

for iSess = 1:length(sessions)
    basename = sessions{iSess};

    cd(fullfile(basepath,basename))
%     
%     sessionInfo = bz_getSessionInfo;
%     spikes = bz_LoadPhy
% bz_GetLFP('all')

%  [clusterIDs, unitQuality, contaminationRate] = maskedClusterQualityKilosort(cd);
%  save([basename '_clusterquality.mat'],'clusterIDs', 'unitQuality','contaminationRate')

% First get info of analogin
load([basename '_analogin'])
% [len_ep, ts_ep, vel_ep, tr_ep, len_ep_fast, ts_ep_fast, vel_ep_fast] =getWheelTrials(analogin)
end

% unc5b env excitation

% edit getWheelTrials