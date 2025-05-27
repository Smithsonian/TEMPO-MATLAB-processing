function success = output_L2_scan_file(fn_in,p_out,expose,NS_bias,fgx,fgy,LAT,LON,H,CITY_VIS,CITY_UV,LIGHTNING,FULL_VIS,FULL_UV,BKGND_VIS,BKGND_UV,VIIRS,...
    bias_vis,bias_uv,bias_lightning,bias_full_vis,bias_full_uv,mu_dx1,mu_dy1,sig_dx1,sig_dy1,success,...
    time,Rsat,Rsun,Rmoon,illum,lunar_zen,lunar_az,sat_zen,sat_az,INDX,STEP,listing)

% Make directory for scan output if it doesn't exist
p_out = fullfile(p_out,'Scan');
if ~exist(p_out,'dir')
    mkdir(p_out);
end

% Name for output L2 file derived from input L1 file name
fn_out = replace(fn_in,'L1','L2');
fn_out = replace(fn_out,'.nc','.mat');
fn_out = replace(fn_out,'G01','');

save(fullfile(p_out,fn_out),'expose','NS_bias','fgx','fgy','LAT','LON','H','CITY_VIS','CITY_UV','LIGHTNING','FULL_VIS','FULL_UV','BKGND_VIS','BKGND_UV','VIIRS',...
    'bias_vis','bias_uv','bias_lightning','bias_full_vis','bias_full_uv','mu_dx1','mu_dy1','sig_dx1','sig_dy1','success',...
    'time','Rsat','Rsun','illum','Rmoon','lunar_zen','lunar_az','sat_zen','sat_az','INDX','STEP','listing')

success = 1;