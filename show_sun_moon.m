function [bore, fh] = show_sun_moon(Rsat,Rsun,Rmoon,lat,lon)

% Sites on Earth
lon = lon*pi/180;
lat = lat*pi/180;

r = 6378.14*[cos(lon(:)).*cos(lat(:)) sin(lon(:)).*cos(lat(:)) sin(lat(:))]';

yaxis = [0; 0; 1];
xaxis = cross(yaxis,Rsat);
xaxis = xaxis/norm(xaxis);

x = zeros(1,size(r,2));
y = x;
for n = 1:size(r,2)
    s = r(:,n) - Rsat;
    s = s/norm(s);
    x(n) = xaxis'*s;
    y(n) = yaxis'*s;
end

fh = figure; plot(asin(x)*180/pi,asin(y)*180/pi,'.')

% Earth disk
re = asin(6378.14/norm(Rsat))*180/pi;
hold on, plot(re*cos((0:360)*pi/180),re*sin((0:360)*pi/180),'k');
axis equal
xlim([-90 90])
ylim([-25 25])

% Boresight (deg)
target_lat = 33.71*pi/180;
target_lon = -96.25267*pi/180; 

r = 6378.14*[cos(target_lon).*cos(target_lat); sin(target_lon).*cos(target_lat); sin(target_lat)];
s = r - Rsat;
s = s/norm(s);
x = xaxis'*s;
y = yaxis'*s;
plot(asin(x)*180/pi,asin(y)*180/pi,'m+')
sb = s;

% Sun
s = Rsun - Rsat;
s = s/norm(s);
x = xaxis'*s;
y = yaxis'*s;
ss = s;

plot(asin(x)*180/pi,asin(y)*180/pi,'ro')

bore = acos(ss'*sb);    % Sun-boresight angle (rad)

% Moon
s = Rmoon - Rsat;
s = s/norm(s);
x = xaxis'*s;
y = yaxis'*s;

plot(asin(x)*180/pi,asin(y)*180/pi,'b*')

xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')

% Solar terminator

r = terminator(Rsun);
m = 0;
cosbeta0 = 6378.14/norm(Rsat);
for n = 1:size(r,2)
    s = r(:,n) - Rsat;
    s = s/norm(s);
    cosbeta = (r(:,n)'*Rsat)/(6378.15*norm(Rsat));
    if (cosbeta>cosbeta0)
        m = m + 1;
        x(m) = xaxis'*s;
        y(m) = yaxis'*s;
    end
end

plot(asin(x)*180/pi,asin(y)*180/pi,'r.')

% Lunar terminator

r = terminator(Rmoon);
m = 0;
x = 0;
y = 0;
for n = 1:size(r,2)
    s = r(:,n) - Rsat;
    s = s/norm(s);
    cosbeta = (r(:,n)'*Rsat)/(6378.15*norm(Rsat));
    if (cosbeta>cosbeta0)
        m = m + 1;
        x(m) = xaxis'*s;
        y(m) = yaxis'*s;
    end
end

x = x(1:m);
y = y(1:m);

plot(asin(x)*180/pi,asin(y)*180/pi,'b.')

% Lunar glint point
lambda0 = 0;
beta = zeros(2,11);
for digit = 1:6
    for n = 1:11
        lambda = lambda0 + (n-1)/10^digit;
        s = lambda*(Rsat-Rmoon) + Rmoon;
        r = 6378.14*s/norm(s);
        beta(1,n)=acosd((Rsat-r)'*r/(6378.14*norm(Rsat-r)));
        beta(2,n)=acosd((Rmoon-r)'*r/(6378.14*norm(Rmoon-r)));
    end
    if (beta(2,1)>beta(1,1))
        k = find(beta(2,:)<=beta(1,:),1);
        lambda0 = lambda0 + (k-2)/10^digit;
    else
        k = find(beta(2,:)>beta(1,:),1);
        lambda0 = lambda0 + (k-2)/10^digit;
    end
end
r = lambda0*(Rsat-Rmoon) + Rmoon;
r = 6378.14*r/norm(r);
s = r - Rsat;
s = s/norm(s);
cosbeta = (r'*Rsat)/(6378.15*norm(Rsat));
if (cosbeta>cosbeta0)
    x = xaxis'*s;
    y = yaxis'*s;
    plot(asin(x)*180/pi,asin(y)*180/pi,'b*')
end
