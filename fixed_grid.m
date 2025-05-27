function [fgx, fgy, r, lv, east] = fixed_grid(lat,lon)

% WGS 84 ellipsoid
WGS84ellipsoid.flat = 1/298.257223563; % Flattening
WGS84ellipsoid.re = 6378.137;          % Earth equatorial radius (km)
WGS84ellipsoid.kappa = (1-WGS84ellipsoid.flat)^2;     % Flattening (derived)
WGS84ellipsoid.rho = 1/(1-WGS84ellipsoid.flat)^2-1;   % Flattening (derived)

% Standard location
a = 42164.17478;
sta = -91*pi/180;

R = a*[cos(sta); sin(sta); 0];

% TEMPO scan orientation
bs_lat = 33.5231*pi/180;
bs_lon = -89.2170*pi/180;

rt = SiteCoordinates(bs_lat,bs_lon,0,WGS84ellipsoid);
st = rt - R;

z = st/norm(st);
x = cross(z,[0; 0; 1]);
x = x/norm(x);
y = cross(z,x);

M = [x y z]';

fgx = zeros(size(lat));
fgy = fgx;

lv = zeros(3,size(lat,1),size(lat,2));
east = lv;
r = lv;

for i = 1:size(lat,1)'
    for j = 1:size(lat,2)

        [r(:,i,j), lv(:,i,j), east(:,i,j)] = SiteCoordinates(lat(i,j)*pi/180,lon(i,j)*pi/180,0,WGS84ellipsoid);

        s = r(:,i,j) - R;
        s = s/norm(s);

        s = M*s;

        fgy(i,j) = -s(2);
        fgx(i,j) = atan2(s(1),s(3));

    end
end