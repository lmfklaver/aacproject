%% These are all the sessions we now used for creating AAC summary plots 
% and AAC Firing Rate plots


dirN  = {...
    'D:/Data/AxoAxo/mouse1/mouse1_180412_2';... %1
    'D:/Data/AxoAxo/mouse1/mouse1_180414';... %2
    'D:/Data/AxoAxo/mouse1/mouse1_180415';... %3
    'D:/Data/AxoAxo/mouse1/mouse1_180501a';...%4
    'D:/Data/AxoAxo/mouse1/mouse1_180501b';...%5
    'D:/Data/AxoAxo/mouse1/mouse1_180502a';...%6
     'D:/Data/AxoAxo/mouse1/mouse1_180502b';...%7
    'D:/Data/AxoAxo/mouse3/mouse3_180627';... %8
    'D:/Data/AxoAxo/mouse3/mouse3_180628';...%9
    'D:/Data/AxoAxo/mouse3/mouse3_180629';...%10
    'D:/Data/AxoAxo/mouse4/mouse4_181114b';...%11
    'D:/Data/AxoAxo/mouse5/mouse5_181112b';...%12
    'D:/Data/AxoAxo/mouse5/mouse5_181116';...%13
    'D:/Data/AxoAxo/mouse6/mouse6_190330';...%14
    'D:/Data/AxoAxo/mouse6/mouse6_190331';...%15
    'D:\Data\Axoaxo\u19_200310_135409';...%16
    'D:\Data\Axoaxo\u19_200313_155505'};...%17
    

% Sess 6 has pulse artifacts
% Sess 7 had the wrong datfile copied over initially and looked empty:
% Still need to check this one!! 
% Sess 10, 11 and 12 have no good epochs %
% Sess 13 has 200ms pulses
% Sess 14 has no AACs


% still check and possibly add:
%     'D:\Data\Axoaxonic_Data_Lianne\u21_200305_153604';...
%     'D:\Data\Axoaxonic_Data_Lianne\u21_200309_142534';...
%     'D:\Data\Axoaxonic_Data_Lianne\u26_200306_172032'};...
%     All the sessions that Kaiser preprocessed under "AAC_Preprocessing_Complete"
%       Earls sessions under "AACRecordingsSorted"

% bad sess:
%     'D:\Data\Axoaxonic_Data_Lianne\u19_200313_120452';...% Pulses are not registered correctly
%     'D:\Data\Axoaxonic_Data_Lianne\m175_200821_151859_2'}  %Arch and ChR both expressed

