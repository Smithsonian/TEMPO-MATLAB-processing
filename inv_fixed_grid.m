function [lat, lon, h] = inv_fixed_grid(fgx,fgy)

% WGS 84 ellipsoid
WGS84ellipsoid.flat = 1/298.257223563; % Flattening
WGS84ellipsoid.re = 6378.137;          % Earth equatorial radius (km)
WGS84ellipsoid.kappa = (1-WGS84ellipsoid.flat)^2;     % Flattening (derived)
WGS84ellipsoid.rho = 1/(1-WGS84ellipsoid.flat)^2-1;   % Flattening (derived)

% Standard location
a = 42164.17478;
sta = -91*pi/180;

fg = WGS84ellipsoid;
fg.a = a;

% R = a*[cos(sta); sin(sta); 0];
R = [a; 0; 0];

% TEMPO scan orientation
bs_lat = 33.5231*pi/180;
bs_lon = -89.2170*pi/180;

rt = SiteCoordinates(bs_lat,bs_lon-sta,0,WGS84ellipsoid);
st = rt - R;

z = st/norm(st);
x = cross(z,[0; 0; 1]);
x = x/norm(x);
y = cross(z,x);

M = [x y z]';

% M0 = [cos(sta) sin(sta) 0; -sin(sta) cos(sta) 0; 0 0 1];


lat = NaN(size(fgx));
lon = fgx;
h = zeros(size(fgx));

for ii = 1:size(fgx,1)
    for jj = 1:size(fgx,2)

        if (isnan(fgx(ii,jj))) continue; end

        s = [ sin(fgx(ii,jj))*sqrt(1-fgy(ii,jj)^2); -fgy(ii,jj); cos(fgx(ii,jj))*sqrt(1-fgy(ii,jj)^2)];

        s = M'*s;

        [~, ~, lat(ii,jj), lon(ii,jj)] = Pierce2(s,fg);

        % Tangent height
        if isnan(lat(ii,jj))

            % Pt. closest to Earth center for first guess
            r = R - (s'*R)*s;
            lat_pt0 = asin(r(3)/norm(r));
            lon_pt0 = atan2(r(2),r(1));
            r_pt0 = SiteCoordinates(lat_pt0,lon_pt0,0,WGS84ellipsoid);
            h_pt0 = norm(r) - norm(r_pt0);

            scale = 1;
            golden_ratio = (1+sqrt(5))/2;
            for nn = 1:100
                d2min = inf;
                d2 = zeros(3,3,3);
                for i = -1:1
                    for j = -1:1
                        for k = -1:1
                            r = SiteCoordinates(lat_pt0+0.01*i*scale,lon_pt0+0.01*j*scale,h_pt0+k*scale,WGS84ellipsoid);
                            v = (r-R) - s'*(r-R)*s;
                            d2(i+2,j+2,k+2) = v'*v;
                            if (d2(i+2,j+2,k+2)<d2min)
                                d2min = d2(i+2,j+2,k+2);
                                i0 = i;
                                j0 = j;
                                k0 = k;
                            end
                        end
                    end
                end
                lat_pt0 = lat_pt0 + 0.01*i0*scale;
                lon_pt0 = lon_pt0 + 0.01*j0*scale;
                h_pt0 = h_pt0 + k0*scale;
                if (i0==0 && j0==0 && k0==0)
                    scale = scale/golden_ratio;
                end
                if (scale<0.001)
                    break
                end
            end

            h(ii,jj) = h_pt0;

        end

    end
end

lon = lon + sta;

lon = lon*180/pi;
lat = lat*180/pi;