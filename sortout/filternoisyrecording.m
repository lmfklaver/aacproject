function buttData = filternoisyrecording(basename)
fileName = [char(basename) '.dat'];
data = bz_LoadBinary(fileName,'frequency',30000,'nChannels',32);
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',30000);
buttData = filtfilt(d,data);
end