function frame = destreak2(frame,qf)

for n = 1:size(frame,2)

    tmp = frame(1:1024,n);
    frame(1:1024,n) = tmp-median(tmp(qf(1:1024,n)==0));
    tmp = frame(1025:2048,n);
    frame(1025:2048,n) = tmp-median(tmp(qf(1025:2048,n)==0));

end