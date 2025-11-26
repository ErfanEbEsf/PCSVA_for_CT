%%
clear all; close all; clc
rng('default')
load ('C:\Users\ebrahimesfahani\Desktop\PCSVA\5_Proposed_PCSVA\SAIST_Denoising_Erfan\Data\TestVal_176_P.mat',...
    'X_test_176')
addpath(genpath('SAIST_Denoising_Erfan'))
%% 
N = size(X_test_176,3);
Sig = 0.05*2;
tic;
for i = 1:4
x = double(X_test_176(:,:,i));
y = x + Sig*randn(size(x));
[xden,Snr,Ssim] = ...
 Image_LASSC_Denoising_Erfan...
 (y,x,Sig);
%X_test_256_den_saist(:,:,i) = xden;
SNR(i) = Snr; SSIM(i) = Ssim;
RMSE(i) = sqrt(mean((xden(:) - x(:)).^2));
end
time = toc
mean_SNR = mean(SNR); std_SNR = std(SNR);
mean_SSIM = mean(SSIM); std_SSIM = std(SSIM);
mean_RMSE = mean(RMSE); std_RMSE = std(RMSE);
SNR = [mean_SNR,std_SNR]
SSIM = [mean_SSIM,std_SSIM]
RMSE = [mean_RMSE,std_RMSE]
%save('Testdata_output_saist_176_2iteration_sig10.mat','SNR','SSIM','RMSE','time')