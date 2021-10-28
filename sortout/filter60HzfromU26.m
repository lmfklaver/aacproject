function out = filter60HzfromU26(basepath)
% 'D:\Data\Axoaxonic_Data_Lianne\u26_200308_144613'
cd(basepath)
basename = bz_BasenameFromBasepath(basepath);
buttData = filternoisyrecording(basename);
out = createDatfromLoadBinary(basename,data,nChans);
end
