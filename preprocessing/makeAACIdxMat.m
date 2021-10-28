

dirN  = {...
    'D:/Data/AxoAxo/mouse1/mouse1_180412_2';... %2
    'D:/Data/AxoAxo/mouse1/mouse1_180414';... %1
    'D:/Data/AxoAxo/mouse1/mouse1_180415';... % 2
    'D:/Data/AxoAxo/mouse1/mouse1_180501a';...%2
    'D:/Data/AxoAxo/mouse1/mouse1_180501b';...%1
    'D:/Data/AxoAxo/mouse1/mouse1_180502a';...%4
    'D:/Data/AxoAxo/mouse1/mouse1_180502b';...%0
    'D:/Data/AxoAxo/mouse3/mouse3_180627';... %0
    'D:/Data/AxoAxo/mouse3/mouse3_180628';...% 0
    'D:/Data/AxoAxo/mouse3/mouse3_180629';...%6
    'D:/Data/AxoAxo/mouse4/mouse4_181114b';...%16
    'D:/Data/AxoAxo/mouse5/mouse5_181112b';...%0
    'D:/Data/AxoAxo/mouse5/mouse5_181116';...%0
    'D:/Data/AxoAxo/mouse6/mouse6_190330';...%0
    'D:/Data/AxoAxo/mouse6/mouse6_190331';...%5
    'D:\Data\Axoaxonic_Data_Lianne\u19_200310_135409';...%3 - but should be 1 AAC, only 2? 6 and 22 definitely not
    %     'D:\Data\Axoaxonic_Data_Lianne\u19_200313_120452';...%2 % should be 4 and 41 - now 4 and 25? Pulses are not registered correctly
    'D:\Data\Axoaxonic_Data_Lianne\u19_200313_155505'};...%3 % correct - 3,13,61
    %     'D:\Data\Axoaxonic_Data_Lianne\u21_200305_153604';...
%     'D:\Data\Axoaxonic_Data_Lianne\u21_200309_142534';...
%     'D:\Data\Axoaxonic_Data_Lianne\u26_200306_172032'};...

% bad sess:
%    
%%
for iSess = [1:6]
    
  cd(dirN{iSess})
  basepath = cd;
  basename = bz_BasenameFromBasepath(cd);
  
[~, ~, aacs] = splitCellTypes(basepath)

aacIDX{iSess} = aacs;
end
save([basename '.modh2l0_5aacs.mat'],'aacIDX')



%12 en 22 niet sess 3
%   13    18    20    21 niet sess 4
% 1     6    17       20    25 niet sess6