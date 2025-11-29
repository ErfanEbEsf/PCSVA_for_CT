clear all, close all, clc
%% Iterative Denoising Rec ADMM
% Maybe MP later
clear all; close all; clc
rng('default')
load ('C:\Users\ebrahimesfahani\Desktop\PCSVA\5_Proposed_PCSVA\SAIST_Denoising_Erfan\Data\TestVal_176_P.mat',...
    'X_test_176')
addpath(genpath('SAIST_Denoising_Erfan')) 
%% Setting parameters 
maxit = 100;
gamma = 25;
rho = 0.05;
alpha = 1;
n = 176;
m = n;
Sig = 0.03;
Sig_added_noise = 0;

%% Defining operators
theta = 0:5*1:179;
Rad =@(x) radon(x,theta);
Rad_adj = @(r) iradon(r,theta(2)-theta(1),n);
%prox2_sigma = @(r,alpha) alpha*r/(alpha+2);
prox_2 = @(u,u0,alpha,rho,gamma)...
u + (1/alpha)*((gamma*rho/(gamma*rho+alpha))-1)*(alpha*Rad_adj(Rad(u))-u0);
%% solution by denoising reconstruction

N = size(X_test_176,3);
X_test_176_rec = zeros(m,n,N);

%XX = X_test_176(:,:,randi([1,N],1,20));
XX = X_test_176;

parfor i = 1:4
x = double(XX(:,:,i));
b = Rad(x+Sig_added_noise*randn(m,n));
u = Rad_adj(b);   
u0 = u;
v = zeros(size(u));
w = zeros(size(u));
SNR_Proposed = zeros(maxit,1);
SSIM_Proposed = zeros(maxit,1);
RMSE_Proposed = zeros(maxit,1);

% Main iterations
for j=1:maxit
    tic;
    u = prox_2(v-(1/rho)*w,u0,alpha,rho,gamma);
    v = Image_LASSC_Denoising_Erfan(u+(1/rho)*w,x,Sig);
    w = w + rho*(u-v);
  
    time(i,j) = toc;
    
    SNR_Proposed(j)=csnr(u,x,0,0);
    SSIM_Proposed(j)=cal_ssim(u,x,0,0);
    RMSE_Proposed(j) = sqrt(mean((u(:) - x(:)).^2));
    fprintf('image= %d iteration= %d SNR= %.2f SSIM= %.4f RMSE=%.4f\n',...
    i,j,SNR_Proposed(j),SSIM_Proposed(j),RMSE_Proposed(j));
 end
 X_test_176_rec(:,:,i) = u;
 Snr(i)=SNR_Proposed(end);Ssim(i)=SSIM_Proposed(end);
 Rmse(i)=RMSE_Proposed(end);
end

mean_SNR = mean(Snr); std_SNR = std(Snr);
mean_SSIM = mean(Ssim); std_SSIM = std(Ssim);
mean_RMSE = mean(Rmse);std_RMSE = std(Rmse);
mean_Time = mean(mean(time));std_Time = std(std(time));
SNR = [mean_SNR,std_SNR]
SSIM = [mean_SSIM,std_SSIM]
RMSE = [mean_RMSE,std_RMSE]
TIME = [mean_Time,std_Time]
% save('Testdata_Recon_Proposed_sr10_2.mat',...
%     'X_test_176_rec','SNR','SSIM','RMSE','time')

% figure;imshow((u));title('Proposed');
% % imwrite(uint8(u),'TGV.png', 'png');
% disparityRange = [0 0.3];
% figure;
% imshow(abs(u-x0),disparityRange);
% title('Error map');
% colormap(gca,jet) 
% colorbar
% 
% figure;imshow((fbp));title('FBP');
% % imwrite(uint8(u),'TGV.png', 'png');
% disparityRange = [0 0.3];
% figure;
% imshow(abs(fbp-x0),disparityRange);
% title('FBP Error map');
% colormap(gca,jet) 
% colorbar