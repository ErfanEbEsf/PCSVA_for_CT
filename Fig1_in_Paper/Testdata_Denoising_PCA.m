%%
clear all; close all; clc
addpath(genpath('C:\Users\ebrahimesfahani\Desktop\8_Examples'))
addpath(genpath(pwd))
rng('default')
load Patient16.mat
%x = X(:,:,170);
Sig = 0.03*1;
x = phantom;
y = x + Sig*randn(256);
y0 = y;
%% 
yy = y;
yy(60:110,100) = 1;
yy(60:110,150) = 1;
yy(60,100:150) = 1;
yy(110,100:150) = 1;
figure, imshow(yy)
 x = (x(60:110,100:150));
 y = (y(60:110,100:150));
%% PCA

[xden_PCA,Snr,Ssim] = ...
 Image_LASSC_Denoising_Erfan...
 (y,x,Sig);
RMSE1 = sqrt(mean((xden_PCA(:) - x(:)).^2));
figure, 
imshow([x, y, xden_PCA])

%% PCA-SVD
%cd('SAIST_Denoising_Erfan_PCA_SVD')
[xden2_PCA_SVD,Snr2,Ssim2] = ...
 Image_LASSC_Denoising_Erfan_PCA_SVD...
 (y,x,Sig);
RMSE2 = sqrt(mean((xden2_PCA_SVD(:) - x(:)).^2));
figure, 
imshow([x, y, xden_PCA, xden2_PCA_SVD])
%% 
figure,
imshow(yy)
figure,
imshow(imresize(y,size(yy)))
%% 
err_PCA = norm(xden_PCA(:) - x(:))/norm(x(:))
err_PCA_SVD = norm(xden2_PCA_SVD(:) - x(:))/norm(x(:))
%% 
noisy_patches = rgb2gray(im2double(imread("noisy_patches.png")));
PCA_patches = rgb2gray(im2double(imread("PCA_patches.png")));
PCA_SVD_patches = rgb2gray(im2double(imread("PCA_SVD_patches.png")));
figure,
imshow([yy,imresize(y,size(yy)),...
  imresize(noisy_patches,size(yy,1)/size(noisy_patches,1)),...
  imresize(PCA_patches,size(yy,1)/size(PCA_patches,1)),...
  imresize(PCA_SVD_patches,size(yy,1)/size(PCA_SVD_patches,1)),...
]')
