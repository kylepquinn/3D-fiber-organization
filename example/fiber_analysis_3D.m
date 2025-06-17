% FIBER_ANALYSIS_3D Perform end-to-end analysis of an intensity image of 
% fibers, including calculation of the pixel-wise fiber orientation,
% segmentation of collagen fibers, and calculating (overall) and (local)
% summary statistics.
%
% Warning: This tutorial is starting to push the limit of requiring a
% relatively large amount of resources and time to process. Keep in mind
% that this tutorial may take some time to run.
%
% See also: calcfibang3D circmean3 directionalVariance3D doublewrite multipageread

clear, clc, close all

%--- Step 1: Load the SHG image ---%
shgImage = multipageread('.\files\fiber_analysis_3D\shg_855.tif');

% Currently, this image is a 16-bit intensity image with a true bit-depth 
% of 13. For most of the analysis in this tutorial, the input image typing 
% (i.e. 8-bit, 16-bit, single, etc) does not matter, except for the
% collagen-positive CNN. For the purpose of this tutorial, the image will
% be converted to have a normalized intensity betweeen 0 and 1.
shgImage = single(shgImage) ./ 8191;


%--- Step 2: Pixel-wise orientation map via `calcfibang3D` ---%
% The function `calcfibang3D` computes the pixel-wise fiber orientation and
% inclination using a vector-weighted summation approach. Along with the 
% input image, the function also expects an XY and Z pixel window size that
% the orientation and inclination will be evaluated over. For example, an 
% XY and Z window value of 7 will calculate the pixel-wise orientation and 
% inclination map over a +/- 7 pixel window about each pixel in the X, Y, 
% and Z directions. Additionally, this function requires a ratio of XY to Z
% resolution to account for a mistmatch between the lateral and axial
% sampling frequencies. The data used in this tutorial has a lateral
% sampling of 1.144 um/pixel, and an axial sampling of 1 um/pixel.
% The output of this function is a pixel-wise orientation angle that ranges
% between [-pi/2 pi/2]. This value is relative to the X/Y axes of the 
% image. For example, an angle of +/- pi/2 means that fiber is pointing 
% parallel to the Y axis, and a value of 0 means that fiber is pointing 
% parallel to the X axis. Negative angles point downward, and positive 
% angles point upward.
% Additionally, the inclination angle is ranged between 0 and pi, where 0
% is pointed along the -z axis, and pi is pointed along the +z axis.
XYtoZRatio = 1.144 / 1;
[orientationImage, inclinationImage] = calcfibang3D(shgImage, 7, 7, XYtoZRatio);

% Flip the inclination. This is not necessary, but easier for me to
% undertand conceptually.
inclinationImage = pi - inclinationImage;

% To display this image, we use an HSV colormap. Generation of a
% pseudo-colored angle map is a multi-step process:
% Firstly, the continuous angle values must be converted to indexed values.
indexedAngles = (orientationImage + pi/2); % Get a positive range of values
indexedAngles = indexedAngles .* (180/pi); % Convert from radians to degrees
indexedAngles = round(indexedAngles) + 1; % Get integer values and add 1 to remove 0 indexing

% Now that the values are indexed values, we can generate the colormap that
% relate to the indexed values.
angleColormap = hsv(181);

% Color the angle image using an indexed colormap
orientationColor = zeros([size(shgImage,[1,2]),3,size(shgImage,3)], 'uint8');
for i = 1:size(shgImage,3)
    tmpColor = ind2rgb(indexedAngles(:,:,i), angleColormap);
    tmpColor = tmpColor .* shgImage(:,:,i);
    orientationColor(:,:,:,i) = uint8(tmpColor.*255); % Convert range from 0-1 to 0-255
end

% Finally, the intensity can be modulated for each pixel based on the input
% image. This also benefits from the image having a normalized intensity
% range between 0 and 1. Please note that is output is **qualitative**. FOr
% future steps, we will directly use the output from `calcfibang2D`.


% The inclination can be processed in a similar fashion, but the
% inclination is usually displayed using the same colormap as the optical
% redox ratio.
inclinationColor = zeros([size(shgImage,[1,2]),3,size(shgImage,3)], 'uint8');
for i = 1:size(shgImage,3)
    inclinationColor(:,:,:,i) = createColorORRMap(inclinationImage(:,:,i), shgImage(:,:,i), 'redoxMin', 0, 'redoxMax', pi, 'redoxInt', pi/180);
end


%--- Step 3: Segment fibers using pretrained collagen-positive CNN ---%
% This is not the only method to segment fibers- other image processing
% techniques can be used in place of this.
%
% See `kylepquinn\collagen-segmentation-cnn` for required resources

% Set up the sigmoid function
sigmoid = @(x) 1./(1+exp(-x));

% Load the pretrained CNN
collagenUNet = importNetworkFromONNX('.\..\unet.onnx');

% Segment each intensity image in the volume
collagenMask = false(size(shgImage));
for i = 1:size(shgImage,3)
    inputImg = dlarray(gpuArray(shgImage(:,:,i)),'SSCB'); % Convert intensity image to array for CNN
    probMap = gather(extractdata(forward(collagenUNet,inputImg))); % Perform forward pass of network with image and get the result
    collagenMask(:,:,i) = (sigmoid(probMap) > 0.56); % Do final softmax and threshold
end


%--- Step 4: Calculate (overall) summary statistics ---%
% Calculate circular summary statstics (mean orientation, mean inclination,
% directional variance) based on each collagen-positive pixel in the image.
[overallMeanOrientation, overallMeanInclination, overallDV] = circmean3(orientationImage, inclinationImage, collagenMask);

% We can also calculate an overall fiber density using the collagen-positive
% mask and a tissue-positive mask. For this tutorial, we are assuming that
% every pixel in the image is tissue-containing.
tissueMask = true(size(collagenMask));
overallFiberDensity = sum(collagenMask(:)) ./ sum(tissueMask(:));


%--- Step 5: Calculate (local) summary statistics ---%
% Calculate the pixelwise circular summary statistics (mean orientation, 
% mean inclination, directional variance) for a local area, determined by a
% box kernel size. The size of the kernel used for this evaulation should 
% be sufficiently large to include multiple fibers within the kernel window.
[localDV, localOrientation, localInclination, localCollagenDensity] = directionalVariance3D(orientationImage, inclinationImage, collagenMask, 40);

% From this local measurement, we can also calculate the local fiber volume
% fraction by using the same function but with a tissue-positive mask
% rather than a collagen-positive mask
[~, ~, ~, localTissueDensity] = directionalVariance3D(0, 0, tissueMask, 40);
localFVF = localCollagenDensity ./ localTissueDensity;


%--- Step 6: Pseudo-color the local summary statistic maps ---%
% For the directional variance and fiber volume fraction maps, we use the
% same colormapping as the redox ratio map.
localDVColor = zeros([size(shgImage,[1,2]),3,size(shgImage,3)], 'uint8');
localFVFColor = zeros([size(shgImage,[1,2]),3,size(shgImage,3)], 'uint8');
for i = 1:size(shgImage,3)
    localDVColor(:,:,:,i) = createColorORRMap(localDV(:,:,i), shgImage(:,:,i).*collagenMask(:,:,i));
    localFVFColor(:,:,:,i) = createColorORRMap(localFVF(:,:,i), shgImage(:,:,i).*collagenMask(:,:,i));
end

% See step 2 for how to pseudo-color orientation maps
indexedAngles = (localOrientation + pi/2); % Get a positive range of values
indexedAngles = indexedAngles .* (180/pi); % Convert from radians to degrees
indexedAngles = round(indexedAngles) + 1; % Get integer values and add 1 to remove 0 indexing

angleColormap = hsv(181);

localOrientationColor = zeros([size(shgImage,[1,2]),3,size(shgImage,3)], 'uint8');
for i = 1:size(shgImage,3)
    tmpColor = ind2rgb(indexedAngles(:,:,i), angleColormap);
    tmpColor = tmpColor .* shgImage(:,:,i) .* collagenMask(:,:,i);
    localOrientationColor(:,:,:,i) = uint8(tmpColor.*255); % Convert range from 0-1 to 0-255
end

localInclinationColor = zeros([size(shgImage,[1,2]),3,size(shgImage,3)], 'uint8');
for i = 1:size(shgImage,3)
    localInclinationColor(:,:,:,i) = createColorORRMap(localInclination(:,:,i), shgImage(:,:,i).*collagenMask(:,:,i), 'redoxMin', 0, 'redoxMax', pi, 'redoxInt', pi/180);
end


%--- Step 7: Quantify local metrics ---%
meanLocalDV = mean(localDV(collagenMask));
sdLocalDV = std(localDV(collagenMask));

meanLocalFVF = mean(localFVF(collagenMask));
sdLocalFVF = std(localFVF(collagenMask));

[localMeanOrientation, localMeanInclination] = circmean3(localOrientation, localInclination, collagenMask);

% Generate a data table with the results.
data = struct('OverallMeanOrientation', overallMeanOrientation * (180/pi),...
              'OverallMeanInclination', overallMeanInclination * (180/pi),...
              'OverallDV', overallDV,...
              'OverallFiberDensity',overallFiberDensity,...
              'LocalMeanOrientation', localMeanOrientation * (180/pi),...
              'LocalMeanInclination', localMeanInclination * (180/pi),...
              'LocalMeanDV', meanLocalDV,...
              'LocalStdDV', sdLocalDV,...
              'LocalMeanFVF', meanLocalFVF,...
              'LocalStdFVF', sdLocalFVF)


%--- Step 8: Save images ---%
% For each pixel-wise map, there should (most of the time) be two output 
% maps: a floating-point image (see `doublewrite` or `hyperwrite`) and a
% pseudo-colored image
%
% For saving 3D stacks, the most reproducible way to write stacks is to
% write the first slice like normal, but then append all other slices. This
% method is overwrite the folder upon writing the first slice, and then add
% on the rest of the slices.
%
% Depending on the speed of the hard drive that these images are being
% saved to, you may encounter this warning:
%
%   Warning: TIFF library error - 'TIFFOpenW:
%   .\files\fiber_analysis_3D\output\inclination.tiff: Cannot
%   open.' - file may be corrupt. 
%
% This warning can be ignored. The images should save properly.

for i = 1:size(shgImage,3)
    if i == 1
        % Orientation
        doublewrite(orientationImage(:,:,i), '.\files\fiber_analysis_3D\output\orientation.tiff');
        imwrite(orientationColor(:,:,:,i), '.\files\fiber_analysis_3D\output\orientation_color.tiff');

        % Inclination
        doublewrite(inclinationImage(:,:,i), '.\files\fiber_analysis_3D\output\inclination.tiff');
        imwrite(inclinationColor(:,:,:,i), '.\files\fiber_analysis_3D\output\inclination_color.tiff');

        % Collagen mask - one of the only exceptions for having one output image.
        % To make an ImageJ-friendly image, the mask should be converted from
        % logical to uint8.
        imwrite(uint8(collagenMask(:,:,i).*255), '.\files\fiber_analysis_3D\output\collagenMask.tiff');
        
        % Local orientation
        doublewrite(localOrientation(:,:,i), '.\files\fiber_analysis_3D\output\localOrientation.tiff');
        imwrite(localOrientationColor(:,:,:,i), '.\files\fiber_analysis_3D\output\localOrientation_color.tiff');

        % Local inclination
        doublewrite(localInclination(:,:,i), '.\files\fiber_analysis_3D\output\localInclination.tiff');
        imwrite(localInclinationColor(:,:,:,i), '.\files\fiber_analysis_3D\output\localInclination_color.tiff');

        % Local directional variance
        doublewrite(localDV(:,:,i), '.\files\fiber_analysis_3D\output\localDV.tiff');
        imwrite(localDVColor(:,:,:,i), '.\files\fiber_analysis_3D\output\localDV_color.tiff');
        
        % Local fiber volume fraction
        doublewrite(localFVF(:,:,i), '.\files\fiber_analysis_3D\output\localFVF.tiff');
        imwrite(localFVFColor(:,:,:,i), '.\files\fiber_analysis_3D\output\localFVF_color.tiff');

    else
        % Orientation
        doublewrite(orientationImage(:,:,i), '.\files\fiber_analysis_3D\output\orientation.tiff','writemode','append');
        imwrite(orientationColor(:,:,:,i), '.\files\fiber_analysis_3D\output\orientation_color.tiff','writemode','append');
        
        % Inclination
        doublewrite(inclinationImage(:,:,i), '.\files\fiber_analysis_3D\output\inclination.tiff','writemode','append');
        imwrite(inclinationColor(:,:,:,i), '.\files\fiber_analysis_3D\output\inclination_color.tiff','writemode','append');

        % Collagen mask - one of the only exceptions for having one output image.
        % To make an ImageJ-friendly image, the mask should be converted from
        % logical to uint8.
        imwrite(uint8(collagenMask(:,:,i).*255), '.\files\fiber_analysis_3D\output\collagenMask.tiff','writemode','append');
        
        % Local orientation
        doublewrite(localOrientation(:,:,i), '.\files\fiber_analysis_3D\output\localOrientation.tiff','writemode','append');
        imwrite(localOrientationColor(:,:,:,i), '.\files\fiber_analysis_3D\output\localOrientation_color.tiff','writemode','append');

        % Local inclination
        doublewrite(localInclination(:,:,i), '.\files\fiber_analysis_3D\output\localInclination.tiff','writemode','append');
        imwrite(localInclinationColor(:,:,:,i), '.\files\fiber_analysis_3D\output\localInclination_color.tiff','writemode','append');

        % Local directional variance
        doublewrite(localDV(:,:,i), '.\files\fiber_analysis_3D\output\localDV.tiff','writemode','append');
        imwrite(localDVColor(:,:,:,i), '.\files\fiber_analysis_3D\output\localDV_color.tiff','writemode','append');
        
        % Local fiber volume fraction
        doublewrite(localFVF(:,:,i), '.\files\fiber_analysis_3D\output\localFVF.tiff','writemode','append');
        imwrite(localFVFColor(:,:,:,i), '.\files\fiber_analysis_3D\output\localFVF_color.tiff','writemode','append');
    end
end