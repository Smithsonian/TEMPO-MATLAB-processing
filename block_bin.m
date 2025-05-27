function hood2 = block_bin(hood,binning)

% Block bin hood
if (binning(1)>1 || binning(2)>1)
    hood2 = zeros(floor(size(hood,1)/binning(2)),floor(size(hood,2)/binning(1)));
    for ii = 1:size(hood2,1)
        for jj = 1:size(hood2,2)
            xhood = hood((ii-1)*binning(2)+(1:binning(2)),(jj-1)*binning(1)+(1:binning(1)));
            hood2(ii,jj) = mean(xhood(:));
        end
    end
else
    hood2 = hood;
end