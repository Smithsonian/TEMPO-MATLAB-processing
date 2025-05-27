% Tile B into A

function [ii, jj, dx, dy, nccx] = align_tile(A,B,max_sites,chip,hood,min_rad)

[~, k] = sort(B(:),'descend');
[i, j] = ind2sub(size(B),k);

% Fill gaps in B
B(~isfinite(B)) = NaN;
B = fill_nan_holes(B')';

num_sites = 0;
m = 0;
d = zeros(1,100);
ii = d;
jj = d;
while (num_sites<max_sites && m<numel(B))

    % Find brightest point in interior neighborhood and <50% overlap

    % Next candidate
    m = m + 1;

    % Test candidate to see if it is interior
    if (i(m)<=(hood(2)/2) || i(m)>(size(B,1)-hood(2)/2) || j(m)<=(hood(1)/2) || j(m)>(size(B,2)-(hood(1)/2)))
        continue
    end

    % Test if too close to others
    d(1:(num_sites-1)) = ((i(m)-ii(1:(num_sites-1)))/3).^2 + (j(m)-jj(1:(num_sites-1))).^2;
    if min(d(1:(num_sites-1)))<(chip(1)/2)^2
        continue
    end

    % Test if bright enough or NaN
    if (B(k(m))<min_rad || isnan(B(k(m))))
        continue
    end

    % Add candidate
    num_sites = num_sites + 1;
    ii(num_sites) = i(m);
    jj(num_sites) = j(m);

end

ii = ii(1:num_sites);
jj = jj(1:num_sites);

% figure, imagesc(B,[0 100]), colormap turbo
% hold on, plot(jj,ii,'r+')

% Loop over sites and match with reference
nccx = zeros(1,num_sites);
dx = nccx;
dy = nccx;
parfor n = 1:num_sites
    template = B(ii(n)+((-chip(2)/2):(chip(2)/2-1)),jj(n)+((-chip(1)/2):(chip(1)/2-1)));
    target = A(ii(n)+((-hood(2)/2):(hood(2)/2-1)),jj(n)+((-hood(1)/2):(hood(1)/2-1)));
    try
        rho = normxcorr2(template,target);
        rho = rho(chip(2):hood(2),chip(1):hood(1));
        [~, k] = max(rho(:));
        [i00, j00] = ind2sub(size(rho),k);
        [ix0, jx0, nccx(n)] = refine3a(i00,j00,rho,1);
        if (abs(ix0-i00)>1 || abs(jx0-j00)>1) nccx(n) = 0; end
        ix = ix0 - (size(rho,1)+1)/2;
        jx = jx0 - (size(rho,2)+1)/2;
        dx(n) = jx;
        dy(n) = -ix;
    catch
        dx(n) = 0;
        dy(n) = 0;
        nccx(n) = 0;
    end
end

