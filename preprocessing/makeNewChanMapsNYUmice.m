%chanmap_thing


chanMap0ind = readNPY('channel_map.npy');
chanMap = chanMap0ind + 1;
positions = readNPY('channel_positions.npy');

kcoords = repmat(1,1,length(chanMap)); % for single shank;
 xcoords = positions(:,1);
 ycoords = positions(:,2);
 connected = kcoords';
 
 save chanMap chanMap chanMap0ind positions kcoords xcoords ycoords connected 
