% MISR + GOES 3D-Winds Code
%
% Original Author: J. Carr, Carr Astronautics, 2018
%                  jcarr@carrastro.com
%
% Methodology: https://www.mdpi.com/2072-4292/10/12/1885
%              https://doi.org/10.3390/rs10121885
%
% Description: Computes site coordinates of a point on earth's surface. 
%
% License: (CC BY 4.0) https://creativecommons.org/licenses/by/4.0/
%
% Work supported through a task on the Support for Atmospheres, Modeling,
% and Data Assimilation (SAMDA) contact at NASA Goddard Space Flight Center
% (NNG17HP01C) awarded to Carr Astronautics, with prime contractor, Science
% Systems and Applications, Inc. (SSAI).

function [r, lv, east, rPhiGC, phiGC] = SiteCoordinates(phiGD,lambda,h,ellipsoid)
%  Compute site coordinates of point on earth's surface.
%  Input:
%    phiGD  == geodetic latitude (rad)
%    lambda == longitude (rad)
%    h      == height above ellipsoid (km)
%    ellipsoid == ellipsoid parameters
%  Output:
%    r=double(3,1) == earth fixed coordinates (km), origin at center of earth, 
%    lv=double(3,1) == local vertical unit vector
%    Z-axis runs to geographic north, and x-axis in the equator runs through prime meridian.
% Earth parameters in ellipsoid structure.
%*******************************************************************

phiGC = atan(ellipsoid.kappa*tan(phiGD));           % Geocentric latitude
rPhiGC = ellipsoid.re/sqrt(1+ellipsoid.rho*sin(phiGC)^2);   % Local radius

coslam = cos(lambda);
sinlam = sin(lambda);
cosphi = cos(phiGC);
sinphi = sin(phiGC);

% Site vector
r = rPhiGC*[cosphi*coslam; cosphi*sinlam; sinphi];

% East
east = [-cosphi*sinlam; cosphi*coslam; 0];

cosphi = cos(phiGD);
sinphi = sin(phiGD);

% Local vertical
lv = [cosphi*coslam; cosphi*sinlam; sinphi];

% Add in height
r = r + h*lv;
