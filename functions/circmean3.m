function varargout = circmean3(varargin)
% CIRCMEAN3 Compute directional summary statstics for 3D spherical data in 
% radians
%
%   [meanTheta, meanPhi] = CIRCMEAN3(thetaVectors, phiVectors) returns the
%   scalar mean direction of orientation in the X-Y plane `meanTheta` and Z
%   inclination `meanPhi` of an array of X-Y plane vectors `thetaVectors`
%   and their associated Z inclination vectors `phiVectors`. The values of 
%   `thetaVectors` are expected to have the range [0 : pi] or 
%   [-pi/2 : pi/2], and `phiVectors` are expected to the the range [0 : pi].
%
%   [meanTheta, meanPhi, directionalVariance] = CIRCMEAN3(...) returns the
%   scalar mean direction of orientation in the X-Y plane `meanTheta` and Z
%   inclination `meanPhi` of an array of X-Y plane vectors `thetaVectors`
%   and their associated Z inclination vectors `phiVectors`, as well as the
%   3D directional variance `directionalVariance`. The value of 
%   `directionalVariance` ranges between 0 (anisotropic) and 1 (isotropic).
%
%   [...] = CIRCMEAN3(thetaVectors, phiVectors, vectorLengths) returns the
%   3D directional summary statistics of an array of X-Y plane vectors 
%   `thetaVectors` and their associated Z inclination vectors `phiVectors`
%   and vector lengths `vectorLengths`.
%
% Editors note (AEW 5/29/25): The convention used in this function is the
% engineering definition of spherical coordinates (i.e. theta === azimuth
% angle, and phi === elevation angle). For the input vector `phiVectors`,
% a value of 0 is orthogonal to the X-Y plane in the +Z direction, a value
% of pi/2 is parallel to the X-Y plane, and a value of pi is orthogonal to
% the X-Y plane is the -Z direction.
%
%
% Written by: Alan E. Woessner (aewoessn@gmail.com)
% Maintained by: Alan E. Woessner (aewoessn@gmail.com) and Kyle P. Quinn
% (kyle.quinn@gmail.com)

% Reshape data
theta = reshape(varargin{1},1,[]);
phi = reshape(varargin{2},1,[]);

if nargin == 2
    mag = 1;
else
    mag = reshape(varargin{3},1,[]);
end

% Calculate C,S,Z - Note we are dealing with bi-directional data so we
% double the angle. 
C = nansum((sin(phi).*cos(2.*theta)).*mag);
S = nansum((sin(phi).*sin(2.*theta)).*mag);
Z = nansum((cos(phi)).*mag);

% Calculate average theta and phi
varargout{1} = atan2(S,C)/2;
varargout{2} = acos(Z/sum(mag));

% Calculate resultant
R = sqrt((C.^2)+(S.^2)+(Z^2));

% Average resultant
rho = R/sum(mag);

% Calculate 3D directional variance
varargout{3} = 1-rho;
end