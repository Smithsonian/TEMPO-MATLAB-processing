% Classify city lights for a full scan after cleaning. Input granules are 
% in folder p_in with file names matching pattern. Output products are put 
% in folder p_out.  The pattern has two parts: {yyyymmdd, Sxxx}
% that identify the the start of scan xxx on the date.  Both are strings.
%
% Status variables returned on output
%
% J. Carr, 2/21/2025

function success = Classify_Scan(p_cleaned, pattern)

success = 0;

% Configuration for processing
configvars_Process_Scan

%% Find cleaned scan belonging to pattern {date, scan}

% Listing will be alphabetized so files are chronologically ordered
listing = dir(fullfile(p_cleaned,'Scan'));
found = 0;
for n = 1:size(listing,1)
    fn = listing(n).name;
    if (contains(fn,pattern{1}) && contains(fn,pattern{2}) && contains(fn,'.mat') && contains(fn,'RADT'))
        load(fullfile(p_cleaned,'Scan',fn),'CITY_VIS')
        found = 1;
        break
    end
end

if (found==0)
    return
end

fn_scan = fn;

% First column in radiance mosaic
i1 = find(isfinite(CITY_VIS(600,:)),1);

%% Available twilight cleaned radiance granules belonging to pattern {date, scan}

% Listing will be alphabetized so files are chronologically ordered
listing_vis = dir(fullfile(p_cleaned,'Granules','VIS'));
listing_uv = dir(fullfile(p_cleaned,'Granules','UV'));

% Find cleaned granule files for granules belonging to scan
belongs = zeros(length(listing_vis),1);
found = 0;
for n = 1:size(listing_vis,1)
    fn = listing_vis(n).name;
    if (contains(fn,pattern{1}) && contains(fn,pattern{2}) && contains(fn,'G01') && contains(fn,'.mat') && contains(fn,'RADT') && ~found)
        % Specified scan beginning on specified date found
        belongs(n) = 1;
        found = 1;      % first granule found
    elseif (contains(fn,pattern{2}) && contains(fn,'.mat') && contains(fn,'RADT') && found)
        % Continues to belong to found scan unless it is granule 1 again
        if (~contains(fn,'G01'))
            belongs(n) = 1;
        else
            % Not part of this instance of the scan number
            break
        end
    elseif (~contains(fn,pattern{2}) && contains(fn,'.mat') && contains(fn,'RADT') && found)
        % Found all granules - done looking
        break
    end
end
if (sum(belongs)==0)
    return
end

%% Prepare spectral library

% Spectral library
warning('off')
lib = readtable(spectral_library);
warning('on')

% Identify categories in library
wvl1 = -inf;
wvl2 = inf;
A = table2array(lib);
lib_wvl = {};
lib_val = {};
lib_cat = {};
lib_labels = lib.Properties.VariableNames;
for j = 3:length(lib_labels)
    cat = strtok(lib_labels(j),'_');
    lib_cat = [lib_cat cat];
    sel = A(:,2)>=wvl_uv(1) & A(:,2)<=wvl_vis(2);
    lib_wvl = [lib_wvl A(sel,2)];
    lib_val = [lib_val A(sel,j)];
    wvl2 = min([wvl2 max(A(sel,2))]);
    wvl1 = max([wvl1 min(A(sel,2))]);
end

% Normalize each library spectrum for unit response = 1 nW/(cm^2 sr) across
% spectral range [wvl1 wvl2] and convert units from photon to energy
for j = 1:length(lib_cat)
    J_per_photon = 6.62607015e-34*299792458./(lib_wvl{j}*1e-9);
    wrad = lib_val{j}.*J_per_photon;
    lambda = lib_wvl{j};
    w = 0;
    for m = 2:length(lambda)
        if(lambda(m)>=wvl1 && lambda(m)<=wvl2)
            w = w + 0.5*(wrad(m)+wrad(m-1))*(lambda(m)-lambda(m-1));
        end
    end
    lib_val{j} = lib_val{j}.*J_per_photon/w;
end

% Unique categories
ulib_cat = unique(lib_cat);

%% Loop over granules belonging to scan and classify bright pixels

% Compositing variables
CONTENTS = NaN(size(CITY_VIS,1),size(CITY_VIS,2),length(ulib_cat));
SIG = NaN(size(CITY_VIS,1),size(CITY_VIS,2));
MU = SIG;

for n = 1:size(listing_vis,1)

    if (~belongs(n))
        continue
    end

    % Read VIS granule variables
    if (use_vis)
        fn_vis = listing_vis(n).name;
        load(fullfile(p_cleaned,'Granules','VIS',fn_vis))
        nframes = size(rad,3);
        rad_vis = rad;
        qf_vis = qf;
        unc_vis = rad_err;
        wvl_vis = wvl;
        clear rad qf rad_err wvl
    else
        rad_vis = [];
        qf_vis = [];
        unc_vis = [];
        wvl_vis = [];
    end
    
    % Read UV granule variables
    if (use_uv)
        fn_uv = listing_uv(n).name;
        load(fullfile(p_cleaned,'Granules','UV',fn_uv))
        nframes = size(rad,3);
        rad_uv = rad;
        qf_uv = qf;
        unc_uv = rad_err;
        wvl_uv = wvl;
        clear rad qf rad_err wvl
    else
        rad_uv = [];
        qf_uv = [];
        unc_uv = [];
        wvl_uv = [];
    end

    % Variables for classifier output
    contents = zeros(size(rad_vis,2),size(rad_vis,3),length(ulib_cat));
    fit_sig2 = zeros(size(rad_vis,2),size(rad_vis,3));
    fit_mu = fit_sig2;
    
    % Loop over frames

    for m = 1:nframes

        % Concatenate bands
        if (use_vis && use_uv)
            wvl = [wvl_uv' wvl_vis'];
            frame = [rad_uv(:,:,m)'/rad_adj(1) rad_vis(:,:,m)'/rad_adj(2)];
            qf = [qf_uv(:,:,m)' qf_vis(:,:,m)'];
            unc = [unc_uv(:,:,m)' unc_vis(:,:,m)'];
        elseif (use_vis)
            wvl = wvl_vis';
            frame = rad_vis(:,:,m)'/rad_adj(2);
            qf = qf_vis(:,:,m)';
            unc = unc_vis(:,:,m)';
        elseif (use_uv)
            wvl = wvl_uv';
            frame = rad_uv(:,:,m)'/rad_adj(1);
            qf = qf_uv(:,:,m)';
            unc = unc_uv(:,:,m)';
        end

        % Truncate to fitting band
        wvl=mean(wvl,1);
        sel = wvl>=wvl1 & wvl<=wvl2;
        wvl = wvl(sel);
        frame = frame(:,sel);
        qf = qf(:,sel);
        unc = unc(:,sel);

        % Make non-nominal radiances NaN
        frame(qf~=0) = NaN;

        % Calculate radiance in VIS band to identify bright pixels to
        % classify
        my_wvl_vis = mean(wvl_vis,2);
        dlambda = mean(my_wvl_vis(2:end)-my_wvl_vis(1:end-1));

        pix_rad = zeros(size(frame,1),1);
        for i = 1:size(pix_rad,1)
            tmp = frame(i,wvl>=wvl_brite(1) & wvl<=wvl_brite(2));
            pix_rad(i) = sum(tmp(~isnan(tmp)))*dlambda;
        end

        % Classify bright pixels
        warning('off')
        for i = 1:size(pix_rad,1)
            if (pix_rad(i)>bright)
                [contents(i,m,:), ~, fit_sig2(i,m), fit_mu(i,m)] = Stepwise_Regression_Classification(wvl,frame(i,:),unc(i,:),lib_cat,lib_wvl,lib_val,ulib_cat,wvl1,wvl2,[]);
            end
        end
        warning('on')

    end 

    % Composite
    i2 = i1 + size(contents,2) - 1;
    CONTENTS(:,i1:i2,:) = contents;
    SIG(:,i1:i2) = sqrt(fit_sig2);
    MU(:,i1:i2) = fit_mu;
    i1 = i2 + 1;

end

% Append classification to the scan file
save(fullfile(p_cleaned,'Scan',fn_scan),'ulib_cat','CONTENTS','SIG','MU','-append')

% Show results
fh1 = figure;
subplot(2,2,1), imagesc(fliplr(squeeze(CONTENTS(:,:,5))),[0 75]), axis off
colormap turbo
colorbar
title(ulib_cat(5))
subplot(2,2,2), imagesc(fliplr(squeeze(CONTENTS(:,:,3))),[0 50]), axis off
colormap turbo
colorbar
title(ulib_cat(3))

subplot(2,2,3), imagesc(fliplr(squeeze(CONTENTS(:,:,2))),[0 30]), axis off
colormap turbo
colorbar
title(ulib_cat(2))
subplot(2,2,4), imagesc(fliplr(squeeze(CONTENTS(:,:,1))),[0 20]), axis off
colormap turbo
colorbar
title(ulib_cat(1))

savefig(fh1,fullfile(p_cleaned,'Figs',[pattern{1} pattern{2}],[pattern{1} pattern{2} '_Classification.fig']))
