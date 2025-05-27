% Process a scan from input folder p_in with file names matching pattern and put
% products in folder p_out.  The pattern has two parts: {yyyymmdd, Sxxx}
% that identify the the start of scan xxx on the date.  Both are strings.
%
% Status variables returned on output
%
% J. Carr, 1/6/2025

function [success, pix_bias, expose, max_col, min_col, err_dx, err_dy, ntp1, illum, avg_lza] = Process_Scan(p_in, p_out, pattern)

success = 0;

% Configuration for processing
configvars_Process_Scan

%% Available twilight radiance files belonging to pattern {date, scan}

% Listing will be alphabetized so files are chronologically ordered
listing = dir(p_in);

% Find RADT netCDF files for granules belonging to scan
belongs = zeros(length(listing),1);
found = 0;
for n = 1:size(listing,1)
    fn = listing(n).name;
    if (contains(fn,pattern{1}) && contains(fn,pattern{2}) && contains(fn,'G01') && contains(fn,'.nc') && contains(fn,'RADT') && ~found)
        % Specified scan beginning on specified date found
        belongs(n) = 1;
        found = 1;      % first granule found
    elseif (contains(fn,pattern{2}) && contains(fn,'.nc') && contains(fn,'RADT') && found)
        % Continues to belong to found scan unless it is granule 1 again
        if (~contains(fn,'G01'))
            belongs(n) = 1;
        else
            % Not part of this instance of the scan number
            break
        end
    elseif (~contains(fn,pattern{2}) && contains(fn,'.nc') && contains(fn,'RADT') && found)
        % Found all granules - done looking
        break
    end
end
if (sum(belongs)==0)
    return
end

%% Process granules

% Compositing variables
CITY_VIS = zeros(2048,canvas_width)-inf;    % City lights in VIS
FULL_VIS = zeros(2048,canvas_width)-inf;    % Full VIS range
BKGND_VIS = zeros(2048,canvas_width)-inf;   % Subtracted VIS background
INDX = zeros(1,canvas_width);               % Index to ID granule
STEP = zeros(1,canvas_width);               % Step within granule
CITY_UV = zeros(2048,canvas_width)-inf;     % City lights in UV
FULL_UV = zeros(2048,canvas_width)-inf;     % Full UV range
BKGND_UV = zeros(2048,canvas_width)-inf;    % Subtracted UV background
LIGHTNING = zeros(2048,canvas_width)-inf;   % Band for lightning
LAT = NaN(2048,canvas_width);
LON = LAT;
TIME = [];
Rsat = [];
Rsun = [];
Rmoon = [];
BKGND_SPEC = [];

first_granule = 1;
max_col = -inf;
min_col = inf;

% Make directory for figures if it doesn't exist
if (~exist(fullfile(p_out,'Figs'),'dir'))
    mkdir(fullfile(p_out,'Figs'));
end
if (~exist(fullfile(p_out,'Figs',[pattern{1} pattern{2}]),'dir'))
    mkdir(fullfile(p_out,'Figs',[pattern{1} pattern{2}]));
end

for n_gran = 1:length(listing)

    if (~belongs(n_gran))
        continue;
    end

    fn = fullfile(p_in,listing(n_gran).name);

    % VIS

    rad = single(ncread(fn,'/band_540_740_nm/radiance'));
    rad_err = single(ncread(fn,'/band_540_740_nm/radiance_error'));
    vis_wvl = ncread(fn,'/band_540_740_nm/nominal_wavelength');
    lat = ncread(fn,'/band_540_740_nm/latitude');
    lon = ncread(fn,'/band_540_740_nm/longitude');
    inr_qf = ncread(fn,'/band_540_740_nm/inr_quality_flag');

    t = ncread(fn,'/time');
    eph_time = ncread(fn,'/inr_input/ephemeris/ephemeris_time');
    sat_X = ncread(fn,'/inr_input/ephemeris/satellite_X');
    sat_Y = ncread(fn,'/inr_input/ephemeris/satellite_Y');
    sat_Z = ncread(fn,'/inr_input/ephemeris/satellite_Z');
    sun_X = ncread(fn,'/geometry/sun_X');
    sun_Y = ncread(fn,'/geometry/sun_Y');
    sun_Z = ncread(fn,'/geometry/sun_Z');
    moon_X = ncread(fn,'/geometry/moon_X');
    moon_Y = ncread(fn,'/geometry/moon_Y');
    moon_Z = ncread(fn,'/geometry/moon_Z');
    rsat = [sat_X sat_Y sat_Z]';
    rsun = [sun_X sun_Y sun_Z]';
    rmoon = [moon_X moon_Y moon_Z]';

    % Interpolate satellite ephemeris to pixel times
    [eph_time, ia] = unique(eph_time);  % protects against repeated samples
    rsat = rsat(:,ia);
    Rsat0 = zeros(3,length(t));
    for n = 1:3
        Rsat0(n,:) = interp1(eph_time',rsat(n,:),t');
    end

    % Gather ephemris information
    TIME = [TIME t'];
    Rsat = [Rsat Rsat0];
    Rsun = [Rsun rsun];
    Rmoon = [Rmoon rmoon];

    % Repair to fix telemtry defects
    [lat, lon] = repair_lat_lon(inr_qf,lat,lon);

    % Clean up
    [rad, qf] = clean_up(rad,do_despeckle(1),do_destreak(1));

    % Mark edge pixels to be neglected
    qf(:,[1:8 2041:end],:) = 4;

    % Background subtraction
    if (do_bkgnd(1)==1)

        bkgnd_rad = zeros(size(rad),'single');

        for j = 1:size(rad,3)

            % Frame with spectrum
            frame = squeeze(rad(:,:,j))';

            % Frame quality flags
            frame_qf = squeeze(qf(:,:,j))';

            % Model background and subtract to leave signal in frame
            [frame, bkgnd_frame] = model_bkgnd4(frame,frame_qf,bkgnd_deg(1),bkgnd_step(1),3);
            for n = 1:refits(1)
                [frame, tmp] = model_bkgnd4(frame,frame_qf,bkgnd_deg(1),bkgnd_step(1),3);
                bkgnd_frame = bkgnd_frame + tmp;
            end

            rad(:,:,j) = frame';
            bkgnd_rad(:,:,j) = bkgnd_frame';

        end

    end

    % Nominal wavelengths
    vis_wvl_nom = single(mean(vis_wvl,2));

    % Spectrum of VIS background
    vis_bkgnd_spec = mean(mean(bkgnd_rad,3),2);

    % Convert radiance from units with ph/s to ones with nW
    J_per_photon = 6.62607015e-34*299792458./(vis_wvl_nom*1e-9);
    for n = 1:size(rad,1)
        rad(n,:,:) = 1e9*rad(n,:,:)*J_per_photon(n);
        bkgnd_rad(n,:,:) = 1e9*bkgnd_rad(n,:,:)*J_per_photon(n);
    end

    % Integrated over spectral ranges to form radiance image
    full = bucket(rad,qf,vis_wvl_nom,[vis_wvl1 vis_wvl2]);
    bkgnd = bucket(bkgnd_rad,zeros(size(qf)),vis_wvl_nom,[vis_wvl1 vis_wvl2]);
    city = bucket(rad,qf,vis_wvl_nom,[city_vis_wvl1 city_vis_wvl2]);

    % Mosaic into scan canvas
    full = fliplr(full);
    city = fliplr(city);
    bkgnd = fliplr(bkgnd);

    lat = fliplr(lat);
    lon = fliplr(lon);

    % Place on canvas comes from scan mechanism position with initial
    % alignemnt from the first pixel of the scan
    if (first_granule==1)
        fn_gran1 = listing(n_gran).name;
        % Read instrument telemetry
        tpix = ncread(fn,'time');
        expose = mean(ncread(fn,'exposure_time'));
        tsma = ncread(fn,'/inr_input/telemetry/mirror/sma_time');
        xsma = ncread(fn,'/inr_input/telemetry/mirror/proc_meas_x');
        ysma = ncread(fn,'/inr_input/telemetry/mirror/proc_meas_y');
        pix_sma = interp1(tsma,xsma,tpix+expose/2)';
        pix_bias = mean(interp1(tsma,ysma,tpix+expose/2));
        % Hack to repair granule
        if isnan(expose)
            expose = 6.2977;
            valid = find(tpix<1e30);
            pix_sma(valid) = interp1(tsma,xsma,tpix(valid)+expose/2)';
            pix_sma = interp1(valid,pix_sma(valid),1:length(pix_sma),'linear','extrap');
            pix_bias = mean(interp1(tsma,ysma,tpix(valid)+expose/2));
            disp(['Repaired: ' listing(nn,:)])
        end
        % Empirically determined coarse alignment to put on canvas
        col_offset = round((3.3973e+04-pix_sma(1))/63.8896)+60+(canvas_width-1230)/2;
        first_granule = 0;
    end

    cols = col_offset + (size(rad,3):-1:1); % canvas columns
    col_offset = col_offset + size(rad,3);

    min_col = min([min_col min(cols)]);
    max_col = max([max_col max(cols)]);

    % Test for off canvas pixels
    if (min(cols)<1 || max(cols)>canvas_width)
        disp(['Warning: pixels off canvas: ' listing(n_gran).name ' cols: [' num2str(min(cols)) ', ' num2str(max(cols)) ']'])
    end
    
    % Fill canvas (supposes it is possible for granules to overlap)  
    for i = 1:2048
        for j = 1:length(cols)
            if (cols(j)<1 || cols(j)>canvas_width)
                continue
            end
            if (city(i,j)>CITY_VIS(i,cols(j)) || isnan(CITY_VIS(i,cols(j))))
                CITY_VIS(i,cols(j)) = city(i,j);
                FULL_VIS(i,cols(j)) = full(i,j);
                BKGND_VIS(i,cols(j)) = bkgnd(i,j);
            end
            INDX(cols(j)) = n_gran;
            STEP(cols(j)) = j;
            LAT(i,cols(j)) = lat(i,j);
            LON(i,cols(j)) = lon(i,j);
        end
    end

    % Save cleaned granule radiance and related data
    output_L1p5_granule_file(listing(n_gran).name,p_out,'VIS',rad,qf,bkgnd,vis_wvl,rad_err);

    disp(['VIS Granule ' listing(n_gran).name])

    % UV

    rad = single(ncread(fn,'/band_290_490_nm/radiance'));
    rad_err = single(ncread(fn,'/band_290_490_nm/radiance_error'));
    uv_wvl = ncread(fn,'/band_290_490_nm/nominal_wavelength');

    % Clean up
    [rad, qf] = clean_up(rad,do_despeckle(2),do_destreak(2));

    % Mark edge pixels to be neglected
    qf(:,[1:8 2041:end],:) = 4;

    % Background subtraction
    if (do_bkgnd(2)==1)

        bkgnd_rad = zeros(size(rad),'single');

        for j = 1:size(rad,3)

            % Frame with spectrum
            frame = squeeze(rad(:,:,j))';

            % Frame quality flags
            frame_qf = squeeze(qf(:,:,j))';

            % Model background and subract to leave signal in frame
            [frame, bkgnd_frame] = model_bkgnd4(frame,frame_qf,bkgnd_deg(2),bkgnd_step(2),3);
            for n = 1:refits(2)
                [frame, tmp] = model_bkgnd4(frame,frame_qf,bkgnd_deg(2),bkgnd_step(2),3);
                bkgnd_frame = bkgnd_frame + tmp;
            end

            rad(:,:,j) = frame';
            bkgnd_rad(:,:,j) = bkgnd_frame';

        end

    end

    % Nominal wavelengths
    uv_wvl_nom = single(mean(uv_wvl,2));

    % Spectrum of UV background
    uv_bkgnd_spec = mean(mean(bkgnd_rad,3),2);

    % Keep background spectrum
    BKGND_SPEC = [BKGND_SPEC; [uv_bkgnd_spec' vis_bkgnd_spec']];

    % Convert radiance from units with ph/s to ones with nW
    J_per_photon = 6.62607015e-34*299792458./(uv_wvl_nom*1e-9);
    for n = 1:size(rad,1)
        rad(n,:,:) = 1e9*rad(n,:,:)*J_per_photon(n);
        bkgnd_rad(n,:,:) = 1e9*bkgnd_rad(n,:,:)*J_per_photon(n);
    end

    % Integrated over spectral ranges to form radiance images
    full = bucket(rad,qf,uv_wvl_nom,[uv_wvl1 uv_wvl2]);
    bkgnd = bucket(bkgnd_rad,zeros(size(qf)),uv_wvl_nom,[uv_wvl1 uv_wvl2]);
    city = bucket(rad,qf,uv_wvl_nom,[city_uv_wvl1 city_uv_wvl2]);

    % Mosaic into scan canvas
    full = fliplr(full);
    city = fliplr(city);
    bkgnd = fliplr(bkgnd);

    for i = 1:2048
        for j = 1:length(cols)
            if (cols(j)<1 || cols(j)>canvas_width)
                continue
            end
            if (city(i,j)>CITY_UV(i,cols(j)) || isnan(CITY_UV(i,cols(j))))
                CITY_UV(i,cols(j)) = city(i,j);
                FULL_UV(i,cols(j)) = full(i,j);
                BKGND_UV(i,cols(j)) = bkgnd(i,j);
            end
        end
    end

    % Integrated over lightning spectral range to form radiance image
    lightning = bucket(rad,qf,uv_wvl_nom,[ln_wvl1 ln_wvl2]);

    % Mosaic lightning into scan canvas
    lightning = fliplr(lightning);

    for i = 1:2048
        for j = 1:length(cols)
            if (cols(j)<1 || cols(j)>canvas_width)
                continue
            end
            if (lightning(i,j)>LIGHTNING(i,cols(j)) || isnan(LIGHTNING(i,cols(j))))
                LIGHTNING(i,cols(j)) = lightning(i,j);
            end
        end
    end

    % Save cleaned granule radiance and related data
    output_L1p5_granule_file(listing(n_gran).name,p_out,'UV',rad,qf,bkgnd,uv_wvl,rad_err);

    disp(['UV Granule ' listing(n_gran).name])

end

% All granules in scan are now cleaned and background subtracted with
% integrated radiances painted on canvas

% Rows can show biases that are calibration and background subtraction
% artifacts that shows as streaks.  This code subtracts them.
bias_vis = NaN;
bias_uv = NaN;
bias_lightning = NaN;
bias_full_vis = NaN;
bias_full_uv = NaN;
if (do_debias(1))
    [CITY_VIS, bias_vis] = debias(CITY_VIS,3);
end
if (do_debias(2))
    [CITY_UV, bias_uv] = debias(CITY_UV,3);
end
if (do_debias(3))
    [LIGHTNING, bias_lightning] = debias(LIGHTNING,3);
end
if (do_debias(4))
    [FULL_VIS, bias_full_vis] = debias(FULL_VIS,3);
end
if (do_debias(5))
    [FULL_UV, bias_full_uv] = debias(FULL_UV,3);
end
if (sum(do_debias)>0)
    fh1 = figure;
    if (do_debias(1))
        hold on, plot(bias_vis,1:2048)
    end
    if (do_debias(2))
        hold on, plot(bias_uv,1:2048)
    end
    if (do_debias(3))
        hold on, plot(bias_lightning,1:2048)
    end
    if (do_debias(4))
        hold on, plot(bias_full_vis,1:2048)
    end
    if (do_debias(5))
        hold on, plot(bias_full_uv,1:2048)
    end
    xlabel('Bias nW/(cm^2 sr)')
    xlim([-5 5])
    ylim([0 2049])
    axis ij
    xlabel('Bias nW/(cm^2 sr)')
    ylabel('Pixels Along Slit')
    title('Row Bias Corrections')
    savefig(fh1,fullfile(p_out,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_Bias.fig']))
end

% Spectrum of subtracted background
fh2 = figure; 
for n = 1:size(BKGND_SPEC,1)
    hold on, plot([uv_wvl_nom' vis_wvl_nom'],BKGND_SPEC(n,:),'.')
end
title('Mean Background Spectra')
xlabel('Wavelength (nm)')
ylabel('Radiance (Ph/(s cm^2 sr nm))')
xlim([290 740])
savefig(fh2,fullfile(p_out,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_BKGND_SPEC.fig']))

% Dark target stats
tmp = CITY_VIS(abs(CITY_VIS(:))<10);
fh3 = figure; histogram(tmp,-10:0.1:10)
mad = median(abs(tmp-median(tmp)));
mu = mean(tmp(abs(tmp-median(tmp))<4*1.4826*mad));
sig = std(tmp(abs(tmp-median(tmp))<4*1.4826*mad));
title(['Dark Background VIS Radiance; Exposure=' num2str(expose) ' s'])
xlabel(['nW/(cm^2 sr), \mu=' num2str(mu,'%4.2f') ', \sigma=' num2str(sig,'%4.2f')])
savefig(fh3,fullfile(p_out,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_Dark_VIS.fig']))
tmp = CITY_UV(abs(CITY_UV(:))<10);
fh4 = figure; histogram(tmp,-10:0.1:10)
mad = median(abs(tmp-median(tmp)));
mu = mean(tmp(abs(tmp-median(tmp))<4*1.4826*mad));
sig = std(tmp(abs(tmp-median(tmp))<4*1.4826*mad));
title(['Dark Background UV Radiance; Exposure=' num2str(expose) ' s'])
xlabel(['nW/(cm^2 sr), \mu=' num2str(mu,'%4.2f') ', \sigma=' num2str(sig,'%4.2f')])
savefig(fh4,fullfile(p_out,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_Dark_UV.fig']))
tmp = LIGHTNING(abs(LIGHTNING(:))<10);
fh5 = figure; histogram(tmp,-10:0.1:10)
mad = median(abs(tmp-median(tmp)));
mu = mean(tmp(abs(tmp-median(tmp))<4*1.4826*mad));
sig = std(tmp(abs(tmp-median(tmp))<4*1.4826*mad));
title(['Dark Background Lightning Radiance; Exposure=' num2str(expose) ' s'])
xlabel(['nW/(cm^2 sr), \mu=' num2str(mu,'%4.2f') ', \sigma=' num2str(sig,'%4.2f')])
savefig(fh5,fullfile(p_out,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_Dark_LIGHTNING.fig']))

% Fix for V3 defect - last navgation of scan may be defective
LAT(:,max_col) = NaN;
LON(:,max_col) = NaN;

% Register to a VIIRS DNB clear-sky mosaic
if (do_registration)

    LAT0 = LAT;
    LON0 = LON;

    % VIIRS Reference Image
    VIIRS = Make_ref_viirs_hires(fn_dnb,LAT0,LON0,[12 6]);
    VIIRS = block_bin(VIIRS,[12 6]);

    % Measure deregistration
    [ii0, jj0, dx0, dy0, nccx0] = align_tile(VIIRS,FULL_VIS,200,chip,hood,min_rad);
    ok0 = nccx0>0.5;    % Tie points to use
    ntp0 = sum(ok0);
    if (ntp0>20)
        % Enough for a bilinear model
        ok0 = madfilt(ok0,dx0,dy0,3);
        mu_dx0 = mean(dx0(ok0));
        mu_dy0 = mean(dy0(ok0));
        sig_dx0 = std(dx0(ok0));
        sig_dy0 = std(dy0(ok0));
        % bilinear warp
        px = bilinear_warp(jj0(ok0),ii0(ok0),dx0(ok0));
        py = bilinear_warp(jj0(ok0),ii0(ok0),dy0(ok0));
    elseif (ntp0>10)
        % Enough for a translation model only
        ok0 = madfilt(ok0,dx0,dy0,3);
        mu_dx0 = mean(dx0(ok0));
        mu_dy0 = mean(dy0(ok0));
        sig_dx0 = std(dx0(ok0));
        sig_dy0 = std(dy0(ok0));
        % translation
        px = [mu_dx0; 0; 0];
        py = [mu_dy0; 0; 0];
    else
        % Insufficient matches for correction
        dx0 = [];
        dy0 = [];
        mu_dx0 = 0;
        mu_dy0 = 0;
        sig_dx0 = NaN;
        sig_dy0 = NaN;
    end

    if (ntp0>10)

        % Reregister
        ii = repmat((1:size(LAT0,1))',1,size(LAT0,2));
        jj = repmat(1:size(LAT0,2),size(LAT0,1),1);
        mdx0 = px(1) + px(2)*jj + px(3)*ii;
        mdy0 = py(1) + py(2)*jj + py(3)*ii;
        LAT = interp2(jj,ii,LAT0,jj+mdx0,ii-mdy0,'linear');
        LON = interp2(jj,ii,LON0,jj+mdx0,ii-mdy0,'linear');

        % Test
        VIIRS = Make_ref_viirs_hires(fn_dnb,LAT,LON,[12 6]);
        VIIRS = block_bin(VIIRS,[12 6]);
        [~, ~, dx1, dy1, nccx1] = align_tile(VIIRS,FULL_VIS,200,chip,hood,min_rad);
        ok1 = nccx1>0.5;
        ntp1 = sum(ok1);
        fh6 = figure;
        if (ntp1>0)

            ok1 = madfilt(ok1,dx1,dy1,3);
            mu_dx1 = mean(dx1(ok1));
            mu_dy1 = mean(dy1(ok1));
            sig_dx1 = std(dx1(ok1));
            sig_dy1 = std(dy1(ok1));
 
            plot(dx0(nccx0>0),dy0(nccx0>0),'r.')
            hold on, plot(dx0(ok0),dy0(ok0),'b.')
            plot(mu_dx0,mu_dy0,'b+')
            plot(mu_dx0+sig_dx0*[-1 1 1 -1 -1],mu_dy0+sig_dy0*[-1 -1 1 1 -1],'b')

            % plot(dx1(nccx1>0),dy1(nccx1>0),'ro')
            plot(dx1(ok1),dy1(ok1),'mo')
            plot(mu_dx1,mu_dy1,'m+')
            plot(mu_dx1+sig_dx1*[-1 1 1 -1 -1],mu_dy1+sig_dy1*[-1 -1 1 1 -1],'m')

            title('Reregistration')
            xlabel('Pixels Cross Slit')
            ylabel('Pixels Along Slit')

            savefig(fh6,fullfile(p_out,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_Reg.fig']))

            if ((abs(mu_dx1)+3*sig_dx1)<1 && (abs(mu_dy1)+3*sig_dy1)<1)
                success = success + 1;
            end

        else
            ntp1 = 0;
            mu_dx1 = 0;
            mu_dy1 = 0;
            sig_dx1 = NaN;
            sig_dy1 = NaN;
        end

    else

        LAT = LAT0;
        LON = LON0;
        ntp1 = 0;
        mu_dx1 = 0;
        mu_dy1 = 0;
        sig_dx1 = NaN;
        sig_dy1 = NaN;

    end

    % Geometric error metrics relative to VIIRS-DNB
    err_dx = abs(mu_dx1)+3*sig_dx1;
    err_dy = abs(mu_dy1)+3*sig_dy1;

    disp(['Done registratiion with VIIRS with error metrics [EW, NS]: [' num2str(err_dx) ', ' num2str(err_dy) '] pixels, N=' num2str(ntp1)])

end

% Navigation variables for TEMPO fixed grid
[fgx, fgy, r, lv, east] = fixed_grid(LAT,LON);
[fgx, fgy] = extend_fg_coords(fgx,fgy);     % extends beyond limb
[~, ~, H] = inv_fixed_grid(fgx,fgy);        % H is height above limb

% Lunar illumination and satellite viewing geometries
lunar_zen = NaN(size(LAT));
sat_zen = lunar_zen;
lunar_az = lunar_zen;
sat_az = lunar_zen;
for i = 1:size(LAT,1)
    for j = min_col:max_col
        V = lv(:,i,j);
        E = east(:,i,j);
        E = E/norm(E);
        N = cross(V,E);
        rmoon = Rmoon(:,j-min_col+1) - r(:,i,j);
        rmoon = rmoon/norm(rmoon);
        lunar_zen(i,j) = acos(rmoon'*V)*180/pi;
        lunar_az(i,j) = atan2(N'*rmoon,E'*rmoon)*180/pi;
        rsat = Rsat(:,j-min_col+1) - r(:,i,j);
        rsat = rsat/norm(rsat);
        sat_zen(i,j) = acos(rsat'*V)*180/pi;
        sat_az(i,j) = atan2(N'*rsat,E'*rsat)*180/pi;
    end
end

% City lights images
fh7 = figure; imagesc(fliplr(CITY_VIS),[0 100]), colorbar, colormap turbo
title(['Radiance Composite: ' pattern{1} ' ' pattern{2} ', band: ' num2str(city_vis_wvl1) ' nm to ' num2str(city_vis_wvl2) ' nm, units: nw/(cm^2 sr)'])
fh8 = figure; imagesc(fliplr(CITY_UV),[0 60]), colorbar, colormap turbo
title(['Radiance Composite: ' pattern{1} ' ' pattern{2} ', band: ' num2str(city_uv_wvl1) ' nm to ' num2str(city_uv_wvl2) ' nm, units: nw/(cm^2 sr)'])
fh9 = figure; imagesc(fliplr(LIGHTNING),[0 60]), colorbar, colormap turbo
title(['Radiance Composite: ' pattern{1} ' ' pattern{2} ', band: ' num2str(ln_wvl1) ' nm to ' num2str(ln_wvl2) ' nm, units: nw/(cm^2 sr)'])
savefig(fh7,fullfile(p_out,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_CITY_VIS.fig']))
savefig(fh8,fullfile(p_out,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_CITY_UV.fig']))
savefig(fh9,fullfile(p_out,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_LIGHTNING.fig']))

% Background images
fh10 = figure; imagesc(fliplr(BKGND_VIS),[0 100]), colorbar, colormap turbo
title(['Background Composite: ' pattern{1} ' ' pattern{2} ', band: ' num2str(vis_wvl1) ' nm to ' num2str(vis_wvl2) ' nm, units: nw/(cm^2 sr)'])
fh11 = figure; imagesc(fliplr(BKGND_UV),[0 300]), colorbar, colormap turbo
title(['Background Composite: ' pattern{1} ' ' pattern{2} ', band: ' num2str(uv_wvl1) ' nm to ' num2str(uv_wvl2) ' nm, units: nw/(cm^2 sr)'])
savefig(fh10,fullfile(p_out,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_BKGND_VIS.fig']))
savefig(fh11,fullfile(p_out,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_BKGND_UV.fig']))

% Sun and moon
try
    elim = isnan(sum(Rsat,1)+sum(Rsun,1)+sum(Rmoon,1));  % eliminate missing epehemeris points
    Rsat(:,elim) = [];
    Rsun(:,elim) = [];
    Rmoon(:,elim) = [];
    mid = round(size(Rmoon,2)/2);
    Rsat0 = Rsat(:,mid);
    Rsun0 = Rsun(:,mid);
    Rmoon0 = Rmoon(:,mid);
    phase = acos(((Rsun0-Rmoon0)'*-Rmoon0)/(norm(Rmoon0)*norm(Rsun0-Rmoon0)))*180/pi;
    illum = 50*(1+cosd(phase));
    [safety, fh12] = show_sun_moon(Rsat0,Rsun0,Rmoon0,LAT([1:16:2048 2048],:),LON([1:16:2048 2048],:));
    title(['SBA = ' num2str(round(safety*180/pi)) ' deg LPA = ' num2str(round(phase)) ' deg Illum. = ' num2str(round(illum)) '%'],'Interpreter','none')
    savefig(fh12,fullfile(p_out,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_Illum.fig']))
catch
    safety = NaN;
    phase = NaN;
    illum = NaN;
    disp(['Unrecoverable anomaly with celestial data: ' pattern{1} pattern{2}])
end

% Average lunar zentih angle for scan
avg_lza = mean(lunar_zen(~isnan(lunar_zen(:))));

success = success + 1;

% Output scan variables
output_L2_scan_file(fn_gran1,p_out,expose,pix_bias,fgx,fgy,LAT,LON,H,CITY_VIS,CITY_UV,LIGHTNING,FULL_VIS,FULL_UV,BKGND_VIS,BKGND_UV,VIIRS,...
    bias_vis,bias_uv,bias_lightning,bias_full_vis,bias_full_uv,mu_dx1,mu_dy1,sig_dx1,sig_dy1,success,...
    TIME,Rsat,Rsun,Rmoon,illum,lunar_zen,lunar_az,sat_zen,sat_az,INDX,STEP,listing);

% Uncomment to close figures inside function
% close([fh1 fh2 fh3 fh4 fh5 fh6 fh7 fh8 fh9 fh10 fh11 fh12])
