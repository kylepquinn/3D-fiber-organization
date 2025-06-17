function varargout = directionalVariance3D(theta,phi,magnitude,kernelSize)
% DIRECTIONALVARIANCE3D Compute local 3D directional variance map
%
%   dv = DIRECTIONALVARIANCE3D(theta, phi, magnitude, kernelSize) 
%   returns the pixel-wise map of local 3D directional variance `dv`, which
%   is calculated based on the pixel-wise map of fiber orientation `theta`
%   (in radians) and fiber inclination `phi` (in radians), and is weighted 
%   based on the map of vector lengths `magnitude`. The local values are 
%   computed based on a specific kernel size `kernelSize`. **Important 
%   note: The `kernelSize` for this function differs compared to 
%   CALCFIBANG3D, in that the local organization map is evaluated over a
%   window of size `kernelSize`.
%
%   [dv, localTheta, localPhi, localDensity] = DIRECTIONALVARIANCE3D(...) 
%   returns the pixel-wise map of local 3D directional variance `dv`, as well
%   as the local orientation maps `localTheta` and `localPhi`, and local 
%   density of collagen-positive pixels `localDensity`.
%
% Editors note (AEW 5/29/25): To allow for rapid evaluation of local 
% directional variance and fiber density, a fast algorithm for summing
% values (see integral images) is used. This drastically speeds up the
% computation, but is strictly limited to square/rectanglular kernel sizes.
%
% See also: calcfibang3D circmean3
%
%
% Written by: Alan E. Woessner (aewoessn@gmail.com)
% Maintained by: Alan E. Woessner (aewoessn@gmail.com) and Kyle P. Quinn
% (kyle.quinn@gmail.com)

% Check for all dimensions of kernel. If they do not exist, then infer
if isscalar(kernelSize)
    kernelSize = [kernelSize,kernelSize,kernelSize];
end

% Convert from circular to cartesian coordinates
xC = single(sin(phi).*cos(2.*theta).*magnitude);
yC = single(sin(phi).*sin(2.*theta).*magnitude);
zC = single(cos(phi).*magnitude);

% Local sums
xCLocal = computeLocalSum3D(xC,kernelSize,'same')./prod(kernelSize); clear xC
yCLocal = computeLocalSum3D(yC,kernelSize,'same')./prod(kernelSize); clear yC
zCLocal = computeLocalSum3D(zC,kernelSize,'same')./prod(kernelSize); clear zC
mLocal = computeLocalSum3D(single(magnitude),kernelSize,'same')./prod(kernelSize);
mLocal(isnan(mLocal)|isinf(mLocal)) = 0;

% Calculate local maps
localTheta = (atan2(yCLocal,xCLocal)./2).*magnitude;
localPhi = acos(zCLocal./mLocal).*magnitude;
R = sqrt((xCLocal.^2)+(yCLocal.^2)+(zCLocal.^2)); clear xCLocal yCLocal zCLocal
dv = 1 - (R ./ (mLocal.*magnitude)); clear R
dv(isnan(dv)|isinf(dv)) = 0;

varargout{1} = dv;
varargout{2} = localTheta;
varargout{3} = localPhi;
varargout{4} = mLocal;
end

function output = computeLocalSum3D(inputImage,kernelSize,shape)
% This function computes a local sum via an integral image

kernelShape = kernelSize;

% Integral image magic. This was pulled from normxcorr2
s = cumsum(padarray(inputImage,kernelShape),1);
c = s(1+kernelShape(1):end-1,:,:)-s(1:end-kernelShape(1)-1,:,:);
s = cumsum(c,2);
c = s(:,1+kernelShape(2):end-1,:)-s(:,1:end-kernelShape(2)-1,:);
s = cumsum(c,3);
output = s(:,:,1+kernelShape(3):end-1)-s(:,:,1:end-kernelShape(3)-1);

if strcmp(shape,'same')
    bound = fix((kernelShape-1)/2);

    if mod(kernelShape(1),2) == 1
        output = output(bound(1)+1:end-bound(1),:,:);
    else
        output = output(bound(1)+2:end-bound(1),:,:);
    end
    
    if mod(kernelShape(2),2) == 1
        output = output(:,bound(2)+1:end-bound(2),:);
    else
        output = output(:,bound(2)+2:end-bound(2),:);
    end

    if mod(kernelShape(3),2) == 1
        output = output(:,:,bound(3)+1:end-bound(3));
    else
        output = output(:,:,bound(3)+2:end-bound(3));
    end
end
end