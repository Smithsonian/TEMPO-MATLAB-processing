function mosaic = fill_nan_holes(mosaic)
% Fill double line holes between tiles where overlap was insufficient
for j = 1:size(mosaic,2)
    for i = 2:size(mosaic,1)-2
        if (isnan(mosaic(i,j))&&isnan(mosaic(i+1,j))&&~isnan(mosaic(i-1,j))&&~isnan(mosaic(i+2,j)))
            mosaic(i,j) = mosaic(i-1,j)+(mosaic(i+2,j)-mosaic(i-1,j))/3;
            mosaic(i+1,j) = mosaic(i-1,j)+2*(mosaic(i+2,j)-mosaic(i-1,j))/3;
        end
    end
end

% Fill single line holes between tiles where overlap was insufficient
for j = 1:size(mosaic,2)
    for i = 2:size(mosaic,1)-1
        if (isnan(mosaic(i,j))&&~isnan(mosaic(i-1,j))&&~isnan(mosaic(i+1,j)))
            mosaic(i,j) = (mosaic(i-1,j) + mosaic(i+1,j))/2;
        end
    end
end

