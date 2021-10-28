basename = 'm175_200821_151859';

for iShank = 1:4
selChannels = chansfromxml(kcoords==iShank)+1;

fnameIn = [basename '.dat'];
fnameOut = [basename '_shank' num2str(shankNum) '.dat'];

% options:
%     - 'duration'  : duration in seconds
%     - 'start'     : start time in seconds
%     - 'channelIx' : subset of channels to be copied


CopyDat(fnameIn,fnameOut,'channelIx', selChannels)
end