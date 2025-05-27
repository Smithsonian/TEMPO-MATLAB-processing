function [frame, sf] = despeckle2(frame,hot,nsig)

sf = zeros(size(frame));

for n = 1:size(frame,1)

    x = abs(frame(n,:)-median(frame(n,:)));
    mad = median(x(~hot(n,:)));

    frame(n,x>nsig*1.4826*mad) = 0;

    sf(n,x>nsig*1.4826*mad) = 2;

end

sf(hot) = 1;