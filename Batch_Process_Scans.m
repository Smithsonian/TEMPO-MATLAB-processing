% Batch process TEMPO nighttime granules
%
% J. Carr, 1/6/2025

p_in = 'F:\TEMPO_DATA_9_months\tempo_data_october_short';
p_out = 'H:\tempo_data_october_short\Cleaned';

start = 20241001;   % yyyymmdd
stop = 20241031;

listing = dir(p_in);
listing = listing(3:end,:);

gran_key = zeros(length(listing),1);
for n = 1:length(listing)
    fn = listing(n).name;
    gran_key(n) = str2num(fn(end-4:end-3));    % Granule number from file name
end

% Open log file
fp = fopen(fullfile(p_out,'processing_log.csv'),'w');
fprintf(fp,'date,scan,success,pix_bias,expose,start_col,stop_col,ew_err,ns_err,ntp,illum,avg_lza,elapsed\n');

for n = 1:length(gran_key)

    % Skip all but granule ones
    if (gran_key(n)~=1)
        continue
    end

    fn = listing(n).name;
    pattern = {fn(19:26),fn(36:39)};    % {yyyymmdd, Sxxx}

    if (str2num(pattern{1})<start || str2num(pattern{1})>stop)
        continue
    end
    
    tic

    % Process scan {yyyymmdd, Sxxx}
    [success, pix_bias, expose, max_col, min_col, err_dx, err_dy, ntp1, illum, avg_lza] = Process_Scan(p_in,p_out,pattern);

    elapsed = toc;

    % Log processing
    fprintf(fp,'%s,%s,%d,%f,%f,%d,%d,%f,%f,%d,%f,%f,%f\n',pattern{1},pattern{2},success,pix_bias,expose,max_col,min_col,err_dx,err_dy,ntp1,illum,avg_lza,elapsed);

    % Close all figures
    close all

    disp(['Elapsed:'  num2str(elapsed) ' s'])

end

fclose(fp);