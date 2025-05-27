% Terminator curve on a spherical Earth illuminated by s

function r = terminator(s)

z = [0; 0; 1];

u = cross(z,s);
u = u/norm(u);
v = cross(s,u);
v = v/norm(v);
w = s/norm(s);

beta = acos(6378.14/norm(s));

phi = (0:360)*pi/180;
r = zeros(3,length(phi));

for n = 1:length(phi)
    r(:,n) = cos(beta)*w + sin(beta)*(cos(phi(n))*u+sin(phi(n))*v);
end

r = r*6378.14;  % km