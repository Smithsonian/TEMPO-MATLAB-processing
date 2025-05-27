function [fgx, fgy] = extend_fg_coords(fgx,fgy)

for m = 1:size(fgx,2)
    sel = ~isnan(fgx(:,m));
    if (sum(sel)<4) continue; end
    t = 1:size(fgx,1);
    p = polyfit(t(sel),fgx(sel,m),3);
    tmp = polyval(p,t);
    fgx(~sel,m) = tmp(~sel);
    p = polyfit(t(sel),fgy(sel,m),3);
    tmp = polyval(p,t);
    fgy(~sel,m) = tmp(~sel);
end