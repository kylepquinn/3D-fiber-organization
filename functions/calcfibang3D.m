function [aSDE,pSDER,bSDER,gSDER] = calcfibang3D(img,kernelSizeXY,kernelSizeZ,XYtoZRatio)
% CALCFIBANG3D Compute pixel-wise 3D orientation map using vector-weighted
% summation algorithm
%
%   [aSDE, pSDER, bSDER, gSDER] = CALCFIBANG3D(img,kernelSizeXY,kernelSizeZ,XYtoZRatio)
%   computes the pixel-wise 3D orientation map of XY orientation `aSDE` and
%   Z inclination `pSDER` of input image volume `img` based on a
%   surrounding window of pixels, dictated by `kernelSizeXY` and
%   `kernelSizeZ`. Additionally, a final (XY resolution / Z resolution)
%   correction factor `XYtoZRatio` is used to account for a mismatch in the
%   sampling frequency in XY relative to Z.
%   **Important note: The pixel size of the window used to calculate the
%   orientation map is evaluated as (2*`kernelSizeXY`)+1 and
%   (2*`kernelSizeZ`)+1.
%
% Editors note (AEW 5/29/25): This function also outputs additional variables `bSDER` and
% `gSDER`, which I am not entirely sure of their function, other than
% computing the 3D directional variance.
%
% See also: calcfibang2D
%
%
% Written by: Zhiyi Liu
% Published in: J Biomed Opt Ex (2015)
% Maintained by: Alan E. Woessner (aewoessn@gmail.com) and Kyle P. Quinn
% (kyle.quinn@gmail.com)

im = img;

im = double(im);
sz = size(im,1);
sz2 = size(im,2);
sz3 = size(im,3);

% Here starts the determination of theta
terimc = zeros(sz,sz2,3*sz3);
terimc(:,:,1:sz3) = im;
terimc(:,:,(sz3+1):(2*sz3)) = im;
terimc(:,:,(2*sz3+1):(3*sz3)) = im;

imc = zeros(sz,sz2,sz3);

for i = 1:sz3
imc(:,:,i) = mean(terimc(:,:,(sz3+i-kernelSizeZ):(sz3+i+kernelSizeZ)),3);
end

clear terimc

meanc = zeros(sz,sz2,sz3);
%h = waitbar(0,'Please Wait');
ij = 0;
nij = sz3;
for i = 1:sz3
    ij = ij+1;
    %waitbar(ij/nij)
    meanc(:,:,i) = calcfibang2D(imc(:,:,i),kernelSizeXY);
end
%close(h)
clear imc

% Here starts the determination of beta
terimj = zeros(3*sz,sz2,sz3);
terimj(1:sz,:,:) = im;
terimj((sz+1):(2*sz),:,:) = im;
terimj((2*sz+1):(3*sz),:,:) = im;

imj = zeros(sz,sz2,sz3);
for i = 1:sz
    imj(i,:,:) = mean(terimj((sz+i-kernelSizeZ):(sz+i+kernelSizeZ),:,:),1);
end
clear terimj
% To calculate the orientation in 2D manner
meanj = zeros(sz,sz2,sz3);
%h = waitbar(0,'Please Wait');
ij = 0;
nij = sz;
for i = 1:sz
    ij = ij+1;
    %waitbar(ij/nij)
    meanj(i,:,:) = calcfibang2D(squeeze(imj(i,:,:)),kernelSizeXY);
end
%close(h)
clear imj
meanj = (pi/2-meanj).*(meanj<=pi/2)+(3*pi/2-meanj).*(meanj>pi/2);

% Here starts the determination of gamma

terimg = zeros(sz,3*sz2,sz3);
terimg(:,1:sz2,:) = im;
terimg(:,(sz2+1):(2*sz2),:) = im;
terimg(:,(2*sz2+1):(3*sz2),:) = im;

img = zeros(sz,sz2,sz3);
for i = 1:sz2
    img(:,i,:) = mean(terimg(:,(sz2+i-kernelSizeZ):(sz2+i+kernelSizeZ),:),2);
end
clear terimg
% To calculate the orientation in 2D manner
meang = zeros(sz,sz2,sz3);
%h = waitbar(0,'Please Wait');
ij = 0;
nij = sz2;
for i = 1:sz2
    ij = ij+1;
    %waitbar(ij/nij)
    meang(:,i,:) = calcfibang2D(squeeze(img(:,i,:)),kernelSizeXY);
end
%close(h)
clear img
meang = (pi/2+meang).*(meang<=pi/2)+(meang-(pi/2)).*(meang>pi/2);

% Here to acquire the orientation of phi
meanp = atan(sqrt(1./(tan(meanj).^2)+1./(tan(meang).^2)));
meanplim = meanp;
for i = 1:sz
    for j = 1:sz2
        for k = 1:sz3
            if meanj(i,j,k) <= pi/2
                meanp(i,j,k) = meanp(i,j,k);
            else
                meanp(i,j,k) = pi-meanp(i,j,k);
            end
        end
    end
end

% Here to acquire the real phi orientation regardless of Z resolution
meanpreal = (pi/2)-(atan(XYtoZRatio*tan((pi/2)-meanplim)));
for i = 1:sz
    for j = 1:sz2
        for k = 1:sz3
            if meanj(i,j,k) <= pi/2
                meanpreal(i,j,k) = meanpreal(i,j,k);
            else
                meanpreal(i,j,k) = pi-meanpreal(i,j,k);
            end
        end
    end
end

% Here to calculate the real beta angle regardless of Z resolution
meanjreal = atan(XYtoZRatio*tan(meanj));
meanjreal = meanjreal.*(meanjreal>=0)+(meanjreal+pi).*(meanjreal<0);

% Here to calculate the real gamma angle regardless of Z resolution
meangreal = atan(XYtoZRatio*tan(meang));
meangreal = meangreal.*(meangreal>=0)+(meangreal+pi).*(meangreal<0);

% Here prepares the final output
aSDE = meanc;
pSDER = meanpreal;
bSDER = meanjreal;
gSDER = meangreal;
end





  
  



