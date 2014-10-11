function [imgRecolor, imgRecolorSim, imgSim] = imgRecolor(imgRGB, type)
%%==========================================================================
% Functionality: image recoloring for enhanced color contrast
%
% This is the example code for the paper
%
% Image Recolorization for the Colorblind
% Jia-Bin Huang, Chu-Song Chen, Tzu-Cheng Jen, and Sheng-Jyh Wang
% Proceedings of IEEE International Conference on Acoustic, Speech, and
% Signal Processing (ICASSP), Taipei, April, 2009
%
% Input:
%   - img: Original RGB image (double)
%   - type: Dicromat type (Protanopia, Deuteranope, Tritanopia)
% Output:
%   - imgSim: Simulated image in RGB
%   - imgRecolor: Recolored image
%   - imgRecolorSim: Simulated recolored image
%
% Jia-Bin Huang
% Electrical and Computer Engineering
% University of Illinois, Urbana-Champaign
% jbhuang1@illinois.edu
% http://www.jiabinhuang.com
%%==========================================================================

%% Testing

[imgHeight, imgWidth, ch] = size(imgRGB);

imgRgbSim = simColorBindImg(imgRGB, type);

%% === Modeling color distribution using GMM ===
% Using smaller size of the image for learning GMM
scale = 0.5;
imgRGBs = imresize(imgRGB, scale, 'bicubic');

% Convert from sRGB to CIELab 
colorTransRgb2Lab = makecform('srgb2lab');
imgLAB = applycform(imgRGBs, colorTransRgb2Lab);

imgL = imgLAB(:,:,1);
imgA = imgLAB(:,:,2);
imgB = imgLAB(:,:,3);

imgLabData = cat(2, imgL(:), imgA(:), imgB(:));

% Fit GMM with color samples with model selection (BIC)
numK = 6;

N = size(imgLabData, 1);
BIC = zeros(1, numK);
gmmObj = cell(1, numK);

for k = 1: numK
    gmmObj{k} = gmdistribution.fit(imgLabData, k, 'CovType', 'diagonal', 'Regularize', 1);
    BIC(k) = gmmObj{k}.BIC;
end

[minBIC, numComponents] = min(BIC);

gmModel = gmmObj{numComponents};

% thetaOrig: parameter of the color distribution of the original image
thetaOrig.mu = gmModel.mu;
thetaOrig.Sigma = gmModel.Sigma;
thetaOrig.weight = gmModel.PComponents;
thetaOrig.K = numComponents;

%% === Compute the simulated color distribution ===

% Assume the covariance matrices of the original and perceived colors are the same
thetaSim = thetaOrig;
thetaSim.mu = simColorLab(thetaOrig.mu, type);

%% === Recolorization by optimizing perceived color contrast ===

% Pre-compute the target contrast
distOrig = pairSymKL(thetaOrig.mu, thetaOrig.Sigma, thetaOrig.K);

% Compute the weight vector
distOrigSim = sqrt(sum((thetaSim.mu - thetaOrig.mu).^2, 2));
weightVec = zeros(1, thetaOrig.K*(thetaOrig.K-1)/2);
ind = 1;
for i = 1: thetaOrig.K -1
    for j = i+1: thetaOrig.K
%         weightVec(ind) = distOrigSim(i) + distOrigSim(j);
        weightVec(ind) = thetaOrig.weight(i) + thetaOrig.weight(j);
        ind = ind + 1;
    end
end

% Initial rotation angles of the key colors
rotInit = atan2(thetaOrig.mu(:, 3), thetaOrig.mu(:, 2));

% Optimization option
opt = optimset( optimset('lsqnonlin') , 'Algorithm', 'levenberg-marquardt', 'MaxIter', 1000,'Diagnostics','off', 'Display','off');

global numComp Sigma origPairDist meanColor typeSim weight;
numComp = thetaOrig.K;
Sigma = thetaOrig.Sigma;
origPairDist = distOrig/sum(distOrig);
meanColor = thetaOrig.mu;
typeSim = type;
weight = weightVec;
% weight = ones(1, thetaOrig.K*(thetaOrig.K-1)/2);

opt = optimset( optimset('lsqnonlin') , 'Algorithm', 'levenberg-marquardt', 'MaxIter', 100,'Diagnostics','off', 'Display','off');

rot = lsqnonlin(@lsq_func, rotInit, [], [], opt);

%% Weighted interpolation using posterior probability

imgLAB = applycform(imgRGB, colorTransRgb2Lab);
imgL = imgLAB(:,:,1); imgA = imgLAB(:,:,2);  imgB = imgLAB(:,:,3);
imgLabData = cat(2, imgL(:), imgA(:), imgB(:));

% Recoloring with linear intepolation
PostProb = posterior(gmModel, imgLabData);
imgTheta = atan2(imgLabData(:,3), imgLabData(:,2));
imgMag = sqrt(imgLabData(:,2).^2 + imgLabData(:,3).^2);
rot = mod(rot, 2*pi);
imgRecolorTheta = imgTheta + PostProb*(rot - rotInit);

imgRecolorLabData = imgLabData;
imgRecolorLabData(:,2) = imgMag.*cos(imgRecolorTheta);
imgRecolorLabData(:,3) = imgMag.*sin(imgRecolorTheta);

imgRecolorL = reshape(imgRecolorLabData(:,1), [imgHeight, imgWidth]);
imgRecolorA = reshape(imgRecolorLabData(:,2), [imgHeight, imgWidth]);
imgRecolorB = reshape(imgRecolorLabData(:,3), [imgHeight, imgWidth]);

imgRecolorLab = cat(3, imgRecolorL, imgRecolorA, imgRecolorB);

% Convert lab to rgb
colorTransLab2Rgb = makecform('lab2srgb');

imgRecolor = applycform(imgRecolorLab, colorTransLab2Rgb);
imgRecolorSim = simColorBindImg(imgRecolor, type);

imgSim = imgRgbSim;

end



function residual = lsq_func(rotTheta)

global numComp Sigma origPairDist meanColor typeSim weight;

magAB = sqrt(meanColor(:,2).^2 + meanColor(:,3).^2);

% Rotation in the A-B plane

newColor = meanColor;
newColor(:,2) = magAB.*cos(rotTheta);
newColor(:,3) = magAB.*sin(rotTheta);

newColorSim = simColorLab(newColor, typeSim);

recolorPairDist  = pairSymKL(newColorSim,  Sigma, numComp);
recolorPairDist = recolorPairDist/sum(recolorPairDist);

residual = weight.*(recolorPairDist - origPairDist);

% disp(['residual = ', num2str(residual)]);

end

function simColor = simColorLab(origColor, type)

colorTransLab2Rgb = makecform('lab2srgb');
colorTransRgb2Lab = makecform('srgb2lab');

% convert key colors to sRGB
imgTempRgb = applycform(origColor, colorTransLab2Rgb);

% simulate perceived colors by the coloblind
imgTempSimRgb = simColorBindImg(imgTempRgb, type);

% convert the simulated colors back to lab

imgTempSimLab = applycform(real(imgTempSimRgb), colorTransRgb2Lab);

simColor = reshape(imgTempSimLab, [size(imgTempSimLab,2) , 3]);

end

function dMat = pairSymKL(mu, Sigma, K)

dMat = zeros(1, K*(K-1)/2);
ind = 1;

for i = 1: K- 1
    for j = i+1: K
        muI = mu(i, :);
        muJ = mu(j, :);
        varI = Sigma(1,:, i);
        varJ = Sigma(1,:, j);
        
        muDiff = muI - muJ;
        varInv = (1./varI) + (1./varJ);
        varTrace = ((varI./varJ) + (varJ./varI)) - 2;
        
        dMat(ind) = sum((muDiff.^2).*varInv + varTrace);
        
        ind = ind + 1;
    end
end

end