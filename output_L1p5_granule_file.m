function success = output_L1p5_granule_file(fn_in,p_out,band,rad,qf,bkgnd,wvl,rad_err)

% Make directory for cleaned granules if it doesn't exist
p_out = fullfile(p_out,'Granules');
if ~exist(p_out,'dir')
    mkdir(p_out);
end
if ~exist(fullfile(p_out,band),'dir')
    mkdir(fullfile(p_out,band));
end

% Name for output L1.5 file derived from input L1 file name
fn_out = replace(fn_in,'L1',['L1p5_' band]);
fn_out = replace(fn_out,'.nc','.mat');

% Background subtracted radinaces could be written but are not to save
% space on disk
% save(fullfile(p_out,fn_out),'rad','qf','bkgnd','wvl','rad_err')

save(fullfile(p_out,band,fn_out),'rad','qf','wvl','rad_err')

success = 1;