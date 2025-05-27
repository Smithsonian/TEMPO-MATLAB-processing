% Stepwise Regression Classifier

function [contents, a, sig2_done, mu_done] = Stepwise_Regression_Classification(lambda,spectrum,uncertainty,lib_cat,lib_wvl,lib_spectrum,ulib_cat,wvl1,wvl2,include)

sel = isnan(spectrum);
lambda(sel) = [];
spectrum(sel) = [];
uncertainty(sel) = [];

sel = lambda>=wvl1 & lambda<=wvl2;
lambda = lambda(sel);
spectrum = spectrum(sel);
uncertainty = uncertainty(sel);

if isempty(lambda)
    contents = NaN(length(ulib_cat),1);
    a = NaN(length(lib_cat),1);
    sig2_done = NaN;
    return
end

A = zeros(length(lambda),size(lib_spectrum,2));

% Interpolate library radiances to
for n = 1:length(lib_spectrum)
    A(:,n) = interp1(lib_wvl{n},lib_spectrum{n},lambda);
end

% Stepwise regression

% in_crowd = [];

for nstep = 1:size(A,2)

    % Test adding new source
    sig2 = inf(1,size(A,2));
    for n = 1:size(A,2)
        if ismember(n,include)
            continue
        end
        [a, r] = wght_spectral_fit(spectrum,uncertainty,A(:,[include n]));
        sig2(n) = r*r'/length(r);
        % Exclude choices with negative light
        if (sum(a<0)>0)
            sig2(n) = inf;
        end
    end

    % Add next most important source
    [minsig, n0] = min(sig2);
    if (~isfinite(minsig))
        % disp(['Stopping at step: ' num2str(nstep)])
        break
    end
    include = [include n0];

end

[a, r] = wght_spectral_fit(spectrum,uncertainty,A(:,include));

sig2_done = r*r'/length(r);
mu_done = mean(r);

% Collect contents by category
contents = zeros(length(ulib_cat),1);
for n = 1:length(include)
    m = find(strcmp(char(lib_cat(include(n))),ulib_cat));
    contents(m) = contents(m) + a(n);
end
