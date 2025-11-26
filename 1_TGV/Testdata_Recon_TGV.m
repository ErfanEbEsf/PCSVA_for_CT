function CT_Rec_TGV
clear all, close all, clc
%% Iterative Denoising Rec Cham_Pock
% Maybe MP later
rng('default')
load ('C:\Users\ebrahimesfahani\Desktop\PCSVA\5_Proposed_PCSVA\SAIST_Denoising_Erfan\Data\TestVal_176_P.mat',...
    'X_test_176')
addpath(genpath('./UTILS'));
%% Setting parameters 
%global alpha0 alpha1 lambda b m n D;
maxit = 500;
sigmaTGV = 1/sqrt(12);
tauTGV = 1/sqrt(12);
lambda = 1e-5;
n = 176;
m = n;
alpha0 = 2;
alpha1 = 0.01;


%% Defining operators
theta = 0:15:179;
Rad =@(x) radon(x,theta);
Rad_adj = @(r) iradon(r,theta(2)-theta(1),n);
prox2_sigma = @(r,lambda) r/(lambda*sigmaTGV+1);
D = @(u) cat(3,dxp(u),dyp(u));
div_1 = @(p) dxm(p(:,:,1)) + dym(p(:,:,2)); 

%% solution by denoising reconstruction
tic,
N = size(X_test_176,3);
X_test_176_rec_TGV = zeros(m,n,N);

XX = X_test_176(:,:,randi([1,N],1,80));

parfor i = 1:N
x = double(X_test_176(:,:,i));
b = Rad(x);

u = Rad_adj(b);        
u_tild = zeros(m,n);
v = D(u);
v_tild = zeros(m,n,2);
p = zeros(m,n,2);
q = zeros(m,n,3);
r = zeros(size(b));
counter=0;

SNR_TGV = zeros(maxit,1);
SSIM_TGV = zeros(maxit,1);
RMSE_TGV = zeros(maxit,1);

% Main iterations
for j=1:maxit
    counter = counter + 1;
    p = projP(p + sigmaTGV*(D(u_tild)-v_tild),alpha1);
    q = projQ(q + sigmaTGV*E(v_tild),alpha0);
    r = prox2_sigma(r + sigmaTGV*(Rad(u_tild) - b),lambda);
    u_old = u;
    u = u + tauTGV*(div_1(p) - Rad_adj(r));
    u = max(0,real(u));
    u_tild = 2*u - u_old;
    v_old = v;
    v = v + tauTGV*(p + div_2(q));
    v_tild = 2*v - v_old;
    
    SNR_TGV(j)=csnr(u,x,0,0);
    SSIM_TGV(j)=cal_ssim(u,x,0,0);
    RMSE_TGV(j) = sqrt(mean((u(:) - x(:)).^2));
    fprintf('image= %d iteration= %d SNR= %.2f SSIM= %.4f RMSE=%.4f\n',...
    i,j,SNR_TGV(j),SSIM_TGV(j),RMSE_TGV(j));
 end
 X_test_176_rec_TGV(:,:,i) = u;
 Snr(i)=SNR_TGV(end);Ssim(i)=SSIM_TGV(end);
 Rmse(i)=RMSE_TGV(end);
end
time = toc
mean_SNR = mean(Snr); std_SNR = std(Snr);
mean_SSIM = mean(Ssim); std_SSIM = std(Ssim);
mean_RMSE = mean(Rmse);std_RMSE = std(Rmse);
SNR = [mean_SNR,std_SNR]
SSIM = [mean_SSIM,std_SSIM]
RMSE = [mean_RMSE,std_RMSE] 
save('Testdata_Recon_TGV_sr15.mat',...
    'X_test_176_rec_TGV','SNR','SSIM','RMSE','time')

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




function z = E(p)
z = zeros(size(p,1),size(p,2),3);
z(:,:,1) = dxm(p(:,:,1));
z(:,:,2) = dym(p(:,:,2));
z(:,:,3) = (dym(p(:,:,1)) + dxm(p(:,:,2)))/2;

function r = div_2(z)
r = zeros(size(z,1),size(z,2),2);
r(:,:,1) = dxp(z(:,:,1)) + dyp(z(:,:,3));
r(:,:,2) = dxp(z(:,:,3)) + dyp(z(:,:,2));


function p = projP(p,alpha1)

  absp = sqrt(abs(p(:,:,1)).^2 + abs(p(:,:,2)).^2);
  denom = max(1,absp/alpha1);
  p(:,:,1) = p(:,:,1)./denom;
  p(:,:,2) = p(:,:,2)./denom;  

function q = projQ(q,alpha0)
  absq = sqrt(abs(q(:,:,1)).^2 + abs(q(:,:,2)).^2 + 2*abs(q(:,:,3)).^2);
  denom = max(1,absq/alpha0);
  q(:,:,1) = q(:,:,1)./denom;
  q(:,:,2) = q(:,:,2)./denom;
  q(:,:,3) = q(:,:,3)./denom;  
 