function lightning = ligntning_detection(rad)

rad(1787,:) = NaN;

x = rad(~isnan(rad(:)) & isfinite(rad(:)));
x = abs(x-median(x));
mad = median(x);
sig = 1.4826*mad;

lightning = rad>6*sig;