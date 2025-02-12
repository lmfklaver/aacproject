%% These are all the sessions we now used for creating AAC summary plots 
% and AAC Firing Rate plots


dirN  = {...
    'D:/Data/Sorting/mouse1/mouse1_180412_2';... %1 ChR2
    'D:/Data/Sorting/mouse1/mouse1_180414';... %2 ChR2
    'D:/Data/Sorting/mouse1/mouse1_180415';... %3 ChR2
    'D:/Data/Sorting/mouse1/mouse1_180501a';...%4 ChR2
    'D:/Data/Sorting/mouse1/mouse1_180501b';...%5 ChR2
    'D:/Data/Sorting/mouse1/mouse1_180502a';...%6 ChR2
    'D:/Data/Sorting/mouse3/mouse3_180627';... %7 ChR2
    'D:/Data/Sorting/mouse3/mouse3_180628';...%8 ChR2
    'D:/Data/Sorting/mouse3/mouse3_180629';...%9 ChR2
    'D:/Data/Sorting/mouse5/mouse5_181112b';...%10 ChR2
    'D:/Data/Sorting/mouse5/mouse5_181116';...%11 ChR2 %200ms stims
    'D:/Data/Sorting/mouse6/mouse6_190331';...%12 ChR2
    'D:\Data\Sorting\u19\u19_200313_155505';...%13 Arch
    'D:/Data\sorting/m217/m217_201027_103818';...%14 Arch
    'D:/Data\sorting/m217/m217_201027_174922';...%15 Arch
    'D:/Data\Sorting/m219/m219_201109_125609';...%16 Arch
    'D:/Data\Sorting/m218/m218_201106_102900';...%17 Arch
    'D:\Data\Sorting\m231\m231\m231_201120_094939';...%18 ChR2 %300ms stims
    'D:\Data\Sorting\m231\m231\m231_201121_130915';...%19 ChR2 %300ms stims
    'D:\Data\Sorting\m418\m418_230523_144934';...%20, CCK Session Arch
    'D:\Data\Sorting\m418\m418_230526_131045';...%21, CCK Session Arch
    'D:\Data\Sorting\m418\m418_230602_140548';...%22, CCK Session Arch
    'D:\Data\Sorting\m420\m420_230609_125733';...%23, CCK Session ChR2
    'F:\AAC_Data\sorting\m305_220421_153424';...% 24, PV Session Arch %560ms
    'F:\AAC_Data\sorting\m306_220502_142219_M';...% 25, PV Session Arch %400ms
    'D:\Data\Sorting\PV3\20160505';...% 26 PV Session ChR2 %200ms, mixed
    'D:\Data\Sorting\PV4\20160225';...% 27 PV Session ChR2 %50ms, mixed
    'D:\Data\Sorting\PV5\20160307';...% 28 PV Session ChR2 %50ms
    'D:\Data\Sorting\PV5\20160308';...% 29 PV Session ChR2 %50ms
    'D:\Data\Sorting\PV5\20160309';...% 30 PV Session ChR2 %50ms
    'D:\Data\Sorting\PV8\20170124';...% 31 PV Session ChR2 %100ms
    'D:\Data\Sorting\PV8\20170125';...% 32 PV Session ChR2 %100ms
    'D:\Data\Sorting\PV8\20170203';...% 33 PV Session ChR2 %500ms mixed
    'D:\Data\Sorting\PV8\20170301'}...% 34 PV Session ChR2 %500ms contaminated
    


    

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

