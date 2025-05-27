% Batch classification of TEMPO nighttime scans
%
% H. Madani, 5/5/2025

p_in = 'H:\tempo_data_october_short\Cleaned';

listing = dir(fullfile(p_in,'Scan','TEMPO*'));
parfor n = 1:size(listing,1)
    fn = listing(n).name;
    pattern = {fn(19:26),fn(36:39)};    % {yyyymmdd, Sxxx}

    Classify_Scan(p_in, pattern)
    elapsed = toc
    elapsed_minutes = elapsed/60
end
