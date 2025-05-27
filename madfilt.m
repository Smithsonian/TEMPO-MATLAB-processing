% Median ABsolute Difference (MAD) filtering

function ok = madfilt(ok,a,b,n)

mada = median(abs(a(ok)-median(a(ok))));
madb = median(abs(b(ok)-median(b(ok))));

ok = ok & abs(a-median(a(ok)))<n*mada*1.4826 & abs(b-median(b(ok)))<n*madb*1.4826;