% Build clearest-sky mosaic
%
% J. Carr, 1/7/2025

% p_in = 'D:\TEMPO Commissioning\Night\Cleaned\Scan';
p_in = 'H:\tempo_data_october_short\Cleaned\Scan';

start = 20241001;   % yyyymmdd
stop = 20241031;

% Configuration for processing
configvars_Process_Scan

% Available scans
listing = dir(p_in);
listing = listing(3:end,:);

% Canvas for mosaic
MOSAIC = NaN(2488,upsample*canvas_width);

%% Loop through scans

for n = 1:length(listing)

    fn = listing(n).name;
    date_str = fn(19:26);

    if (str2num(date_str)<start || str2num(date_str)>stop)
        continue
    end

    load(fullfile(p_in,fn),'success','FULL_VIS','LIGHTNING','fgx','fgy','lunar_zen')

    % Skip scans that were not completely cleaned and registered to VIIRS
    if (success<2)
        continue
    end

    scan = FULL_VIS;    % Thing to mosaic

    % Mask out moonlit pixels
    if (moonless)
        scan(lunar_zen<90) = -inf;
    end

    % Floor for city lights radiance to avoid stacking noise
    scan(scan<5) = 0;

    % Mask lightning
    lightning = ligntning_detection(LIGHTNING);
    scan(lightning) = -inf;

    % Indices for the fixed-grid places of pixels
    ii = round(-fgy/fg0_dy + floor(2048/2+220)+0.5);
    jj = round(fgx/(fg0_dx/upsample) + floor(upsample*canvas_width/2)+0.5);

    for m = 1:numel(ii)
        if (ii(m)>0 && jj(m)>0 && ii(m)<=2488 && jj(m)<=upsample*canvas_width)
            MOSAIC(ii(m),jj(m)) = max([scan(m) MOSAIC(ii(m),jj(m))]);
        end
    end

end

figure, imagesc(MOSAIC,[0 100]), colormap turbo, colorbar, title('TEMPO Radiance (nW/(cm^2 sr nm))')

% Standard fixed grid definition
fgx = repmat((fg0_dx/upsample)*((1:upsample*canvas_width)-(floor(upsample*canvas_width/2)+0.5)),2488,1);
fgy = repmat(fg0_dy*(floor(2048/2+220)+0.5-(1:2488)'),1,upsample*canvas_width);
[fg_LAT, fg_LON, H] = inv_fixed_grid(fgx,fgy);

% VIIRS Reference Image
VIIRS = Make_ref_viirs_hires(fn_dnb,fg_LAT,fg_LON,[12/upsample 6]);
VIIRS = block_bin(VIIRS,[12/upsample 6]);

figure, imagesc(VIIRS,[0 100]), colormap turbo, colorbar, title('VIIRS-DNB Radiance (nW/(cm^2 sr nm))')

% Registration check
[~, ~, dx1, dy1, nccx1] = align_tile(VIIRS,MOSAIC,200,chip,hood,min_rad);
ok1 = nccx1>0.5;
ntp1 = sum(ok1);

ok1 = madfilt(ok1,dx1,dy1,3);
mu_dx1 = mean(dx1(ok1));
mu_dy1 = mean(dy1(ok1));
sig_dx1 = std(dx1(ok1));
sig_dy1 = std(dy1(ok1));

figure, plot(dx1(ok1),dy1(ok1),'mo')
hold on, plot(mu_dx1,mu_dy1,'m+')
plot(mu_dx1+sig_dx1*[-1 1 1 -1 -1],mu_dy1+sig_dy1*[-1 -1 1 1 -1],'m')
title('Registration')
xlabel('Pixels Cross Slit')
ylabel('Pixels Along Slit')

% Radiance comparison
viirs=block_bin(VIIRS,[4 4]);
tempo=block_bin(MOSAIC,[4 4]);
figure, density_histogram(viirs(:),tempo(:),[5 2 97],[5 2 97],5:4:97,5:4:97,1,1)
xlabel('VIIRS-DNB (nW/(cm^2 sr nm)')
ylabel('TEMPO (nW/(cm^2 sr nm)')
title('Radiance (nW/(cm^2 sr nm))')

sel = ~isnan(tempo(:))&viirs(:)>5&tempo(:)>5&viirs(:)<100&tempo(:)<100&abs(viirs(:)-tempo(:))<20;
figure, histogram(viirs(sel)./tempo(sel))
title('Ratio of Collocated Similar Pixels')
xlabel('VIIRS/TEMPO')

save('Clearest_Sky_Mosaic_oct_1scan.mat','MOSAIC','fg_LAT','fg_LON')