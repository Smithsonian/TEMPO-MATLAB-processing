% MISR + GOES 3D-Winds Code
%
% Original Author: J. Carr, Carr Astronautics, 2018
%                  jcarr@carrastro.com
%
% Methodology: https://www.mdpi.com/2072-4292/10/12/1885
%              https://doi.org/10.3390/rs10121885
%
% Description: Performs interpolation on the correlation surface.
%
% License: (CC BY 4.0) https://creativecommons.org/licenses/by/4.0/
%
% Work supported through a task on the Support for Atmospheres, Modeling,
% and Data Assimilation (SAMDA) contact at NASA Goddard Space Flight Center
% (NNG17HP01C) awarded to Carr Astronautics, with prime contractor, Science
% Systems and Applications, Inc. (SSAI).

%%%
%
% Quadratic surface fit with uniqueness metric.
%
% Copyright (C) 2013 Carr Astronautics Co.
%
% Classification:  None
%
% Author: J. Carr, 10/23/2013
%
%%%

function [i0, j0, rho_max, v1, v2, uniqueness] = refine3a(i0,j0,RHO,quicker)

v1 = zeros(2,1); v2 = v1;
uniqueness = 0;

% shift of of peak
row_shift = i0 - floor(size(RHO,1)/2) - 1;
col_shift = j0 - floor(size(RHO,2)/2) - 1;

% Test if maximum was found at an interior address in correlation matrix
if (i0>1 && i0<size(RHO,1) && j0>1 && j0<size(RHO,2))
    
    %     % FOR ACCELERATED NCC ONLY ==>
    %     % Added to account for the fact that accelerated NCC really outputs sgn(RHO)*RHO^2
    %     R = sign(RHO((i0-1):(i0+1),(j0-1):(j0+1)))...
    %         .*sqrt(abs(RHO((i0-1):(i0+1),(j0-1):(j0+1))));
    R = RHO((i0-1):(i0+1),(j0-1):(j0+1));
    
    % Maximum is in interior -- interpolation is possible.
    % Find biquadratic model for correlation matrix in neighborhood of maximum using
    % least-squares fit of max. and its 8 nearest neighbors.
    % Model is rho(x,y)=A0+Ax*x+Ay*y+(Axx*x*x+2*Axy*x*y+Ayy*y*y)/2
    
    X = [-1 -1 -1; 0 0 0; 1 1 1]; Y = X';
    
    A0 = (5/9)*trace(R'*ones(3,3))-(1/3)*trace(R'*(X.*X))-(1/3)*trace(R'*(Y.*Y));
    Ax = (1/6)*trace(R'*X);
    Ay = (1/6)*trace(R'*Y);
    Axy = (1/4)*trace(R'*(X.*Y));
    Axx = -(2/3)*trace(R'*ones(3,3)) + trace(R'*(X.*X));
    Ayy = -(2/3)*trace(R'*ones(3,3)) + trace(R'*(Y.*Y));
    
    % Determinant of second derivative matrix {{Axx,Axy},{Axy,Ayy}}
    DetA = Axx*Ayy - Axy*Axy;
    
    % Test for saddle surface
    if (DetA > 0)
        % Surface is not a saddle surface.
        % Find correction to integer address of correlation maximum
        ddx = (Ay*Axy - Ax*Ayy)/DetA;
        ddy = (Ax*Axy - Ay*Axx)/DetA;
        %Add correction to correlation maximum location
        row_shift = row_shift + ddx;
        col_shift = col_shift + ddy;
        % Correlation value is that of surface maximum
        rho_max = A0 - 0.5*(Axx*ddx*ddx+2*Axy*ddx*ddy+Ayy*ddy*ddy);
        % Spectral decomposition of 2nd derivative matrix gives error ellipse
        radical = sqrt((Axx-Ayy)*(Axx-Ayy)+4*Axy*Axy);
        % Find eigenvalues -- optimized to minimize round-off error
        eigenvalue1 = 0.5*(Axx + Ayy - radical);
        eigenvalue2 = DetA/eigenvalue1;
        % Find eigenvectors
        if (Axx < Ayy)
            v1(1) = Ayy-Axx+radical;
            v1(2) = -2*Axy;
            v2(1) = -v1(2);
            v2(2) = v1(1);
        elseif (Axx ~= Ayy || Axy ~= 0.0)
            v2(1) = Axx-Ayy+radical;
            v2(2) = 2*Axy;
            v1(1) = -v2(2);
            v1(2) = v2(1);
        else
            % Special case for multiple of identity matrix
            v1(1) = 1;
            v1(2) = 0;
            v2(1) = 0;
            v2(2) = 1;
        end
        % Scale eigenvectors inversely to square-root of -eigenvalues.
        % Eigenvectors form error ellipse semi-axes.
        % Note that eigenvalues are all negative.
        scale = 1/sqrt(-eigenvalue1*(v1(1)*v1(1)+v1(2)*v1(2)));
        v1 = v1*scale;
        scale = 1.0/sqrt(-eigenvalue2*(v2(1)*v2(1)+v2(2)*v2(2)));
        v2 = v2*scale;
        % else
        %   % Interpolation correction is too large.
        % 	% No confidence in result => set correlation to zero.
        % 	rho_max = 0.5;
        % end
    else
        % Peak resembles a saddle surface when fit, cannot interpolate
        rho_max = 0;
    end
    
    if (~quicker)
        
        % Measure uniqueness of peak w.r.t. other high values
        peak = RHO(i0,j0);
        
        % Slower uniqueness algorithm
        second = min(RHO(:));
        for i = 2:size(RHO,1)-1
            for j = 2:size(RHO,2)-1
                if (i==i0 && j==j0) continue; end
                local_max = ...
                    RHO(i,j) > max([RHO(i-1,j) RHO(i+1,j) RHO(i,j-1) RHO(i,j+1) ...
                    RHO(i-1,j-1) RHO(i-1,j+1) RHO(i+1,j-1) RHO(i+1,j+1)]);
                if (local_max)
                    if (RHO(i,j) > second)
                        second = RHO(i,j);
                    end
                end
            end
        end
        uniqueness = peak/second;
        
    end
    
    % Refined position of peak
    i0 = row_shift + floor(size(RHO,1)/2) + 1;
    j0 = col_shift + floor(size(RHO,2)/2) + 1;
    
else
    
    % Maximum is on boundary -- interpolation is impossible.
    % No confidence in result => set correlation to zero.
    rho_max = 0;
    
end
