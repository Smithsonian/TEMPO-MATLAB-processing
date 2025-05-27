% Locate point where line-of-sight pierces the ellipsoid

function [r, lv, Lat, Lon, OnEarth, ECA] = Pierce2(los,fg)

re = fg.re;
kappa = fg.kappa;
a = fg.a;

% Determine if Earth is pierced by LOS
q1 = los(1)*los(1) + los(2)*los(2) + los(3)*los(3)/kappa;
q2 = los(1)*a/re;
q3 = (a*a)/(re*re) - 1.0;
discrim = q2*q2 - q1*q3;
if (discrim > 0.0)
    OnEarth = 1;
    % Calculate slant range to the site
    slantrange = -re*(q2 + sqrt(discrim))/q1;
    % Compute site position vector
    r1 = a + slantrange*los(1);
    r2 = slantrange*los(2);
    r3 = slantrange*los(3);
    r = [r1; r2; r3];
    ECA = acos(r1/norm(r));
    % Convert to geodetic latitude and longitude relative to station
    Lat = atan(r3/(kappa*sqrt(r1*r1+r2*r2)));
    Lon = atan(r2/r1);
    % Local vertical at site
    lv = [cos(Lat)*cos(Lon); cos(Lat)*sin(Lon); sin(Lat)];
else
    r = NaN(3,1);
    lv = r;
    Lat = NaN;
    Lon = NaN; 
    OnEarth = 0;
    ECA = NaN;
end
