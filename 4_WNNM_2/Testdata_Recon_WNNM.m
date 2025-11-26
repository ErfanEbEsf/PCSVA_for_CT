clear all, close all, clc
%% Iterative Denoising Rec Cham_Pock
% Maybe MP later
clear all; close all; clc
rng('default')
load ('C:\Users\ebrahimesfahani\Desktop\PCSVA\5_Proposed_PCSVA\SAIST_Denoising_Erfan\Data\TestVal_176_P.mat',...
    'X_test_176')
addpath(genpath('SAIST_Denoising_Erfan')) 
%% Setting parameters 
maxit = 100;
sigma = 1/10;
tau = 1;
lambda = 1e-5;
n = 512;
m = n;
Sig = 0.03;
Sig_added_noise = 0;

%% Defining operators
theta = 0:5:179;
Rad =@(x) radon(x,theta);
Rad_adj = @(r) iradon(r,theta(2)-theta(1),n);
prox2_sigma = @(r,lambda) r/(lambda*sigma+1);

%% solution by denoising reconstruction

N = size(X_test_176,3);
X_test_176_rec = zeros(m,n,N);

XX = X_test_176(:,:,randi([1,N],1,20));
%XX = X_test_176;

parfor i = 1:20
x = double(XX(:,:,i));
x = imresize(x,[n,n]);
b = Rad(x);
u = Rad_adj(b);   
u_tild = zeros(m,n);
v_tild = zeros(m,n,2);
r = zeros(size(b));
SNR_Proposed = zeros(maxit,1);
SSIM_Proposed = zeros(maxit,1);
RMSE_Proposed = zeros(maxit,1);

% Main iterations
for j=1:maxit
    r = prox2_sigma (r + sigma*(Rad(u_tild) -b),lambda);
    u_old = u;
    u = u + tau*(- Rad_adj(r));   
    tic;
    u = Image_LASSC_Denoising_Erfan(u,x,Sig);
    Time(i,j) = toc;
    u_tild = 2*u - u_old;
    
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
mean_TIME = mean(mean(Time));std_TIME=std(std(Time));

SNR = [mean_SNR,std_SNR]
SSIM = [mean_SSIM,std_SSIM]
RMSE = [mean_RMSE,std_RMSE] 
TIME = [mean_TIME,std_TIME]
%save('Testdata_Recon_WNNM_sr10.mat',...
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