function [rad, bias] = debias(rad,nsig)

bias = zeros(size(rad,1),1);

for i = 1:size(rad,1)

    line = rad(i,:);

    test = line(line>-inf); 
    city = zeros(size(test));
    z = abs(test-median(test));
    mad = median(z);
    city(z>nsig*1.4826*mad) = 1;

    bias(i) = mean(test(~city));

    rad(i,:) = rad(i,:) - bias(i);

end