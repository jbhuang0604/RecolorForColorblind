function imgSim = simColorBindImg(img, type)

%%==========================================================================
% Functionality: simulated images perceived by the colorblind.
%
% The simulation algorithm is described in
% H. Brettel, F. Viénot, and J. D. Mollon (1997)
% Computerized Simulation of Color Appearance for Dichromats. J. Opt. Soc. Am. A 14, 2647 - 2655.
%
% Input:
%   - img: Original RGB image (double)
%   - type: Dicromat type (Protanopia, Deuteranope, Tritanopia)
% Output:
%   - imgSim: Simulated image in RGB
%
% Jia-Bin Huang
% Electrical and Computer Engineering
% University of Illinois, Urbana-Champaign
% jbhuang1@illinois.edu
% http://www.jiabinhuang.com
%%==========================================================================

sizeImg = size(img); 
if(length(sizeImg)==3)
    imgHeight = sizeImg(1);
    imgWidth  = sizeImg(2);
    imgR = img(:,:,1);
    imgG = img(:,:,2);
    imgB = img(:,:,3);
else
    imgHeight = 1;
    imgWidth  = sizeImg(1);
    imgR = img(:,1);
    imgG = img(:,2);
    imgB = img(:,3);
end
GAMMA  = 2.2;

imgRGBVec = cat(2, imgR(:), imgG(:), imgB(:))';

% Remove sRGB gamma nonlinearity
imgRGBVec = imgRGBVec.^GAMMA;

% Transformation matrix RGB -> LMS
rgb2lms = [17.8824 43.5161 4.11935; 3.45565 27.1554 3.86714; 0.0299566 0.184309 1.46709] ;
% Transformation matrix LMS -> RGB
lms2rgb = [0.0809 -0.1305 0.1167; -0.0102 0.0540 -0.1136; -0.0004 -0.0041 0.6935];

% Map rgb image to lms color space
imgLMSVec = rgb2lms*imgRGBVec;

switch type
    case 'Protanopia'
        % Projection matrix for Protanopia in the LMS space
        T = [0 2.02344 -2.52581; 0 1 0; 0 0 1] ;
    case 'Deuteranope'
        T = [1 0 0; 0.494207 0 1.24827; 0 0 1] ;
    case 'Tritanopia'
        T = [1 0 0; 0 1 0; -0.395913 0.801109 0] ;
end

% Compute LMS values perceived by the colorblind
imgSimLMS = T* imgLMSVec;

% Convert from LMS to RGB
imgSimRGBVec = lms2rgb*imgSimLMS;

imgSimR = imgSimRGBVec(1,:);
imgSimG = imgSimRGBVec(2,:);
imgSimB = imgSimRGBVec(3,:);

imgSimR = imgSimR.^(1/GAMMA);
imgSimG = imgSimG.^(1/GAMMA);
imgSimB = imgSimB.^(1/GAMMA);

imgSimR = reshape(imgSimR, [imgHeight, imgWidth]);
imgSimG = reshape(imgSimG, [imgHeight, imgWidth]);
imgSimB = reshape(imgSimB, [imgHeight, imgWidth]);

imgSim = cat(3, imgSimR, imgSimG, imgSimB);

imgSim = real(imgSim);