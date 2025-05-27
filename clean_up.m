function [rad, qf, bad, hot_pix] = clean_up(rad,do_speckle,do_destreak)

% Quality flags
qf = zeros(size(rad,1),size(rad,2),size(rad,3),'int8');

% Find anomalously oversubtracted quadrants
mquad = zeros(2,size(rad,3));
for n = 1:size(rad,3)
    frame = squeeze(rad(:,:,n))';
    tmp = frame(1:1024,:);
    mquad(1,n) = mean(tmp(:));
    tmp = frame(1025:2048,:);
    mquad(2,n) = mean(tmp(:));
end
x = abs(mquad-median(mquad,2));
mad = median(x,2);
bad = abs(x)>6*1.4826*mad;

% Find persistently hot pixels
hot_pix = zeros(size(rad,2),size(rad,1));
for n = 1:size(rad,3)
    frame = squeeze(rad(:,:,n))';
    [~, sf] = despeckle(frame,10);
    hot_pix = hot_pix + sf;
end

% Despeckle/Destreak frames
parfor n = 1:size(rad,3)
    frame0 = squeeze(rad(:,:,n))';
    if (do_speckle)
        [~, qf0] = despeckle2(frame0,hot_pix>3,6);
    else
        qf0 = hot_pix>3;
    end
    if (do_destreak)
        frame = destreak2(frame0,qf0);
    else
        frame = frame0;
    end
    qf(:,:,n) = qf0';
    frame(qf0~=0) = 0;
    rad(:,:,n) = frame';
end

% Apply qf to bad quadrants
for n = 1:size(rad,3)
    if (bad(1,n))
        qf(:,1:1024,n) = 3;
    end
    if (bad(2,n))
        qf(:,1025:2048,n) = 3;
    end
end

% if (~do_bkgnd)
%     return
% end

% % Model and subtract background
% for n = 1:size(rad,3)
%     frame0 = squeeze(rad(:,:,n))';
% for j = 1:size(frame0,2)
%     [city(1:1024,j), bkgnd(1:1024,j)] = model_bkgnd(frame0(1:1024,j),2);
%     [city(1025:2048,j), bkgnd(1025:2048,j)] = model_bkgnd(city(1025:2048,j),2);
% end
% end