%%
clear all; close all; clc
rng('default')
load ('C:\Users\ebrahimesfahani\Desktop\PCSVA\5_Proposed_PCSVA\SAIST_Denoising_Erfan\Data\AAPM_test.mat')
addpath(genpath('UTILS')) 
%% 
Op = 0;
FS = 40;
N = size(X,3);
nx = 512;
tic;
for i = 1:N
x = double(X(:,:,i));
x = imresize(x,[nx,nx]);
y = double(Y(:,:,i));
y = imresize(y,[nx,nx]);
SSIM_Noisy(i) = cal_ssim(y,x,0,0);
RMSE_Noisy(i) = sqrt(mean((y(:) - x(:)).^2));

%Sig = 0.03; %sqrt(iedd(y))
[xden] = TGV_den(y);
%X_test_256_den_saist(:,:,i) = xden;
SNR(i) = csnr(x,xden); SSIM(i) = cal_ssim(x,xden);
RMSE(i) = sqrt(mean((xden(:) - x(:)).^2));


x = rgb2gray(insertText(x,[1,1],['Ground Truth',newline,'SSIM= 100%',newline, 'RMSE= 0%'],'BoxOpacity',Op,'FontSize',FS,TextColor='white'));
y = rgb2gray(insertText(y,[1,1],['Low-dose Noisy',newline, 'SSIM= ',num2str(100*SSIM_Noisy(i),'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_Noisy(i),'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));
xden = rgb2gray(insertText(xden,[1,1],['Proposed',newline, 'SSIM= ',num2str(100*SSIM(i),'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE(i),'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));
z=[y,xden,x];
figure,imshow(z)
% imwrite(z,['AAPM',num2str(i),'.png'])
end
    time = toc    
mean_SNR = mean(SNR); std_SNR = std(SNR);
mean_SSIM = mean(SSIM); std_SSIM = std(SSIM);
mean_RMSE = mean(RMSE); std_RMSE = std(RMSE);
SNR = [mean_SNR,std_SNR]
SSIM = [mean_SSIM,std_SSIM]
RMSE = [mean_RMSE,std_RMSE]

mean_SSIM_Noisy = mean(SSIM_Noisy); std_SSIM_Noisy = std(SSIM_Noisy);
mean_RMSE_Noisy = mean(RMSE_Noisy); std_RMSE_Noisy = std(RMSE_Noisy);
SSIM_Noisy = [mean_SSIM_Noisy,std_SSIM_Noisy]
RMSE_Noisy = [mean_RMSE_Noisy,std_RMSE_Noisy]
%save('Testdata_output_saist_176_2iteration_sig10.mat','SNR','SSIM','RMSE','time')