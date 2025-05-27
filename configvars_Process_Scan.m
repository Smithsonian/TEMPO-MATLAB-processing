% Folder for radiance files

% Clean up flags [vis uv]
do_despeckle = [1 1];
do_destreak = [0 0];
do_bkgnd = [1 1];
bkgnd_deg = [3 3];
bkgnd_step = [1 1];
refits = [1 1];
do_debias = [1 1 1 1 1];    % [city_vis city_uv lightning full_vis full_uv]

do_registration = 1;

% VIIRS-DNB registration
fn_dnb = 'SVDNB_npp_20221201-20221231_75N180W_vcmslcfg_v10_c202301111000.avg_rade9h.tif';
chip = [32 96];             % Template [EW NS] size in pixels
hood = chip + [20 40];      % Search area for matching
min_rad = 20;   % nW/(cm^2 sr)  % Minimum template radiance

vis_wvl1 = 540; % VIS full range
vis_wvl2 = 740; % nm

city_vis_wvl1 = 540; % City lights VIS range
city_vis_wvl2 = 640;

uv_wvl1 = 290;  % UV full range
uv_wvl2 = 490;

city_uv_wvl1 = 390; % City lights UV range
city_uv_wvl2 = 490;

ln_wvl1 = 310;  % Lightning range
ln_wvl2 = 340;

canvas_width = 1630; %1230; 

% Standard fixed grid sampling distance in rad
fg0_dx = 123e-6;
fg0_dy = 41.49e-6;
upsample = 2;

% Compositing with moonlight
moonless = 1;

% City lights classifier
use_vis = 1;
use_uv = 1;
wvl_vis = [540 740]; % nm
wvl_uv = [290 490];  % nm
wvl_brite = [540 640];  % nm
bright = 5;     % nW/(cm^2 sr)
rad_adj = [1.10 1.06]; % [uv vis]
spectral_library = 'euro_noaa_lights_database_unscaled.csv';