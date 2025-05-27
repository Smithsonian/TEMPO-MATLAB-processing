function [lat, lon] = repair_lat_lon(inr_qf,lat,lon)

if (sum(inr_qf(:)==2)==0)
    return
end

for i = 1:size(lat,1)
    x = 1:size(lat,2);
    y1 = lat(i,:);
    y2 = lon(i,:);
    remove = inr_qf(i,:)==2;
    x(remove) = [];
    y1(remove) = [];
    y2(remove) = [];
    lat(i,:) = interp1(x,y1,1:size(lat,2),'linear','extrap');
    lon(i,:) = interp1(x,y2,1:size(lon,2),'linear','extrap');
end
