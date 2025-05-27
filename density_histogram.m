function density_histogram(x,y,xrng,yrng,xticks,yticks,do_log,do_diag)

js = floor((x-xrng(1))/xrng(2))+1;
is = floor((y-yrng(1))/yrng(2))+1;

cnt = zeros(floor((yrng(3)-yrng(1))/yrng(2)),floor((xrng(3)-xrng(1))/xrng(2)));
for n = 1:length(x)
    if (x(n)<xrng(1) || x(n)>=xrng(3) || y(n)<yrng(1) || y(n)>=yrng(3))
        continue
    end
    try
        cnt(is(n),js(n)) = cnt(is(n),js(n)) + 1;
    catch
        continue
    end
end

if (do_log)
    cnt(cnt<=0) = NaN;
    cnt = log10(cnt);
    max_cnt = max(cnt(:));
    imagesc(cnt,[0 max_cnt]), axis xy, axis square, colormap turbo
    colorbar('Ticks',0:floor(max_cnt),'TickLabels',10.^(0:floor(max_cnt)))
else
    imagesc(cnt), axis xy, axis square, colormap turbo, colorbar
end

ax = gca;
ax.XAxis.TickValues = floor((xticks-xrng(1))/xrng(2))+0.5;
for n = 1:length(ax.XAxis.TickValues)
    ax.XAxis.TickLabels{n} = num2str(xticks(n));
%     ax.XAxis.TickLabels{n} = num2str((str2num(ax.XAxis.TickLabels{n})-1)*xrng(2)+xrng(1));
end
ax.YAxis.TickValues = floor((yticks-yrng(1))/yrng(2))+0.5;
for n = 1:length(ax.YAxis.TickValues)
    ax.YAxis.TickLabels{n} = num2str(yticks(n));
%     ax.YAxis.TickLabels{n} = num2str((str2num(ax.YAxis.TickLabels{n})-1)*yrng(2)+yrng(1));
end

if (do_diag)
    hold on, plot([0.5 size(cnt,2)+0.5],[0.5 size(cnt,1)+0.5],'w')
end