function [signal, bkgnd, p] = model_bkgnd4(frame,dqf,degree,step,nsig)

% Eliminate NaNs
tmp = frame;
tmp(isnan(tmp)) = 0;

% Integrate over spectral bins 
test = sum(tmp,2);

% Find bright locations (cities) in top and bottom of frame w/ nsig MAD
if (step)
    city_top = zeros(1024,1);
    z = abs(test(1:1024)-median(test(1:1024)));
    mad = median(z);
    city_top(z>nsig*1.4826*mad) = 1;
    city_bot = zeros(1024,1);
    z = abs(test(1025:2048)-median(test(1025:2048)));
    mad = median(z);
    city_bot(z>nsig*1.4826*mad) = 1;
    city = [city_top; city_bot];
else
    city = zeros(2048,1);
    z = abs(test(1:2048)-median(test(1:2048)));
    mad = median(z);
    city(z>nsig*1.4826*mad) = 1;
end

% Independent variable for fit
x = ((1:size(tmp,1))-size(tmp,1)/2+0.5)/size(tmp,1);

% Allocate memory for background frame
bkgnd = NaN(size(tmp));

% Model coefficients
p = zeros(degree+1+step,size(tmp,2));

% Design matrix for fitting with polynomial and optionally a step
F = zeros(length(x),degree+1+step);
for n = 0:degree
    F(:,n+1) = x.^n;
end
if (step)
    F(1:1024,degree+2) = -0.5;
    F(1025:2048,degree+2) = 0.5;
end

% Background fit for each spectral bin
for j = 1:size(tmp,2)

    % Use only nominal data quality and not bright locations
    sel = ~city & dqf(:,j)==0;

    % Skip background fitting if not enough useful data
    if (sum(sel)<0.5*length(sel))
        continue
    end

    % Fit
    A = F(sel,:)'*F(sel,:);
    v = F(sel,:)'*tmp(sel,j);
    p(:,j) = A\v;

    % Evaluate fit for the background
    bkgnd(:,j) = F*p(:,j);

end

% Subtract background from frame
signal = frame - bkgnd;