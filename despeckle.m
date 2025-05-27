function [frame, sf] = despeckle(frame,nsig)

sf = zeros(size(frame));

for n = 1:size(frame,1)

    x = abs(frame(n,:)-median(frame(n,:)));
    mad = median(x);

    frame(n,abs(x)>nsig*1.4826*mad) = 0;

    sf(n,abs(x)>nsig*1.4826*mad) = 1;

end