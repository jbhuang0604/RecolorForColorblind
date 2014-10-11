%% Demo
% 
% This is the reference implementation of the paper
%
% Image Recolorization for the Colorblind
% Jia-Bin Huang, Chu-Song Chen, Tzu-Cheng Jen, and Sheng-Jyh Wang
% Proceedings of IEEE International Conference on Acoustic, Speech, and
% Signal Processing (ICASSP), Taipei, April, 2009
%
% If you find the software useful, please consider cite our paper.
% 
% Input:
%   - img: Original RGB image (double)
%   - type: Dicromat type (Protanopia, Deuteranope, Tritanopia)
% Output:
%   - imgSim: Simulated image in RGB
%   - imgRecolor: Recolored image
%   - imgRecolorSim: Simulated recolored image

% Jia-Bin Huang
% Electrical and Computer Engineering
% University of Illinois, Urbana-Champaign
% jbhuang1@illinois.edu
% http://www.jiabinhuang.com

%% 
imagePath = 'images';
resultPath = 'results';

imageName = 'flower.jpg';

type = 'Protanopia';
imgRGB = im2double(imread(fullfile(imagePath, imageName)));

% Three types of color vision deficiency.
%      'Protanopia'
%      'Deuteranope'
%      'Tritanopia'

tic;
[imgRecolorRgb, imgRecolorSimRgb, imgSimRgb] = imgRecolor(img, type);
t = toc;
disp(['Image recoloring done in ', num2str(t), ' seconds']);

% Write images
imwrite(imgRecolorRgb, fullfile(resultPath, [imageName(1:end-4), '_Recolor_',type, '.png']));
imwrite(imgRecolorSimRgb, fullfile(resultPath, [imageName(1:end-4), '_RecolorSim_', type,'.png']));
imwrite(imgSimRgb, fullfile(resultPath, [imageName(1:end-4), '_OrigSim_', type, '.png']));
imwrite(img, fullfile(resultPath, [imageName(1:end-4), '_Orig.png']));

%% Show images 
figure(1)
subplot(2,2,1), imshow(img);
title('Original image');
subplot(2,2,2), imshow(imgSimRgb);
title(['Simulation of the origianl image: (', type,')']);
subplot(2,2,3), imshow(imgRecolorRgb);
title('Recolored image');
subplot(2,2,4), imshow(imgRecolorSimRgb);
title(['Simulation of the recolored image: (', type,')']);
