%%
clear all; close all; clc
rng('default')
load ('D:\content from internal HDD\Erfan\codes\Machine Learnng Datasets\Ferdowsi\FUMPE\Images_mfile\Processed Data\TrainTestVal_176_P.mat',...
    'X_test_176')
addpath(genpath('./UTILS'));
addpath(genpath('D:\content from internal HDD\Erfan\codes\BM3D_MRI_toolbox\Utils\BM3D'))
%% 
N = size(X_test_176,3);
Sig = 0.1;
tic;
parfor i = 1:N
x = double(X_test_176(:,:,i));
y = x + Sig*randn(size(x));
[xden] = TGV_den(y);
%X_test_256_den_saist(:,:,i) = xden;
SNR(i) = csnr(x,xden); SSIM(i) = cal_ssim(x,xden);
RMSE(i) = sqrt(mean((xden(:) - x(:)).^2));
end
time = toc
mean_SNR = mean(SNR); std_SNR = std(SNR);
mean_SSIM = mean(SSIM); std_SSIM = std(SSIM);
mean_RMSE = mean(RMSE); std_RMSE = std(RMSE);
SNR = [mean_SNR,std_SNR]
SSIM = [mean_SSIM,std_SSIM]
RSME = [mean_RMSE,std_RMSE]
% save('Testdata_output_TGV_176_sig10.mat','SNR','SSIM','RMSE','time')