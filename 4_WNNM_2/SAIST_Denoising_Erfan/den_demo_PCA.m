clc, close all, clear all
addpath(genpath('./Data'))
rng('default')
Sigma=15/255;
    x=double(imread('barbara.tif'))/255;
    y=x+randn(size(x))*Sigma;
    x_saist=Image_LASSC_Denoising_Erfan(y,x,sqrt(iedd2(y)));
figure,
subplot(1,2,1),imshow((y))
subplot(1,2,2),imshow((x_saist))