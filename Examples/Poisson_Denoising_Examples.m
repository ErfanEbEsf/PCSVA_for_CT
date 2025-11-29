clear all, close all, clc
load ('C:\Users\ebrahimesfahani\Desktop\PCSVA\5_Proposed_PCSVA\SAIST_Denoising_Erfan\Data\AAPM_test.mat')
rng('default')
%figure, imshow([x3, Seg])
% figure, imshow([imresize(x3,176/512),imresize(Seg,176/512)])
%% Noisy
n = 512; % size(x,1);
m = n;
i=3;
x = double(X(:,:,i));
x = imresize(x,[n ,n ]);
y = double(Y(:,:,i));
 
%% Proposed
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\5_Proposed_PCSVA')
addpath(genpath('SAIST_Denoising_Erfan'))
Sig = 0.035;
tic;
[Prop,Snr,Ssim] = ...
 Image_LASSC_Denoising_Erfan...
 (y,x,Sig);
SNR_Prop = Snr; SSIM_Prop = Ssim;
RMSE_Prop = sqrt(mean((Prop(:) - x(:)).^2));
time = toc
Sig = 0.03;
%% WNNM
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\4_WNNM_2')
addpath(genpath('SAIST_Denoising_Erfan'))
tic;
[WNNM,Snr,Ssim] = ...
 Image_LASSC_Denoising_Erfan...
 (y,x,Sig);
SNR_WNNM = Snr; SSIM_WNNM = Ssim;
RMSE_WNNM = sqrt(mean((WNNM(:) - x(:)).^2));
time = toc

%% SAIST
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\3_SAIST')
addpath(genpath('SAIST_Denoising_Erfan'))
tic;
[SAIST,Snr,Ssim] = ...
 Image_LASSC_Denoising_Erfan...
 (y,x,Sig);
SNR_SAIST = Snr; SSIM_SAIST = Ssim;
RMSE_SAIST = sqrt(mean((SAIST(:) - x(:)).^2));
time = toc

%% BM3D
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\2_BM3D')
addpath(genpath('.\Utils'))
tic;
[~,u] = BM3D2(1,y,255*Sig);
u = double(u);
SNR_BM3D=csnr(u,x,0,0);
SSIM_BM3d= cal_ssim(u,x,0,0);
RMSE_BM3d= sqrt(mean((u(:) - x(:)).^2));
time = toc
BM3d = u;
%% TGV
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\1_TGV')
addpath(genpath('./UTILS'));
tic;
[TGV] = TGV_den(y);
SNR_TGV = csnr(x,TGV); SSIM_TGV = cal_ssim(x,TGV);
RMSE_TGV = sqrt(mean((TGV(:) - x(:)).^2));
time = toc

%% Classical PCA
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\6_classical_PCA')
addpath(genpath('SAIST_Denoising_Erfan'))
tic;
[PCA,Snr,Ssim] = ...
 Image_LASSC_Denoising_Erfan...
 (y,x,Sig);
SNR_PCA = Snr; SSIM_PCA = Ssim;
RMSE_PCA = sqrt(mean((PCA(:) - x(:)).^2));
time = toc

%% Output
SSIM_Noisy = cal_ssim(y,x,0,0);
RMSE_Noisy = sqrt(mean((y(:) - x(:)).^2));
Delta = 100;
x1 = 230; x2 = x1 + Delta; %70 %120
y1 = 250; y2 = y1 + Delta; %78 %117 
Op = 0;
FS = 40;

BB = zeros(1,size(y,2));
z = [y;BB;TGV;BB;PCA;BB;BM3d;BB;...
    SAIST;BB;WNNM;BB;Prop;BB;x];
d = abs(z - [x;BB;x;BB;x;BB;x;BB;...
    x;BB;x;BB;x;BB;x]);


x = rgb2gray(insertText(x,[1,1],['Ground Truth',newline,'SSIM= 100%',newline, 'RMSE= 0%'],'BoxOpacity',Op,'FontSize',FS,TextColor='white'));
y = rgb2gray(insertText(y,[1,1],['Noisy',newline, 'SSIM= ',num2str(100*SSIM_Noisy,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_Noisy,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));
PCA = rgb2gray(insertText(PCA,[1,1],['PCA',newline, 'SSIM= ',num2str(100*SSIM_PCA,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_PCA,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));
TGV = rgb2gray(insertText(TGV,[1,1],['TGV',newline, 'SSIM= ',num2str(100*SSIM_TGV,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_TGV,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white')); 
SAIST = rgb2gray(insertText(SAIST,[1,1],['SAIST',newline, 'SSIM= ',num2str(100*SSIM_SAIST,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_SAIST,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));
BM3d = rgb2gray(insertText(BM3d,[1,1],['BM3D',newline, 'SSIM= ',num2str(100*SSIM_BM3d,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_BM3d,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));
WNNM = rgb2gray(insertText(WNNM,[1,1],['WNNM',newline, 'SSIM= ',num2str(100*SSIM_WNNM,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_WNNM,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));
Prop = rgb2gray(insertText(Prop,[1,1],['Proposed',newline, 'SSIM= ',num2str(100*SSIM_Prop,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_Prop,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));

zz = [y;BB;TGV;BB;PCA;BB;BM3d;BB;...
    SAIST;BB;WNNM;BB;Prop;BB;x];

BB = zeros(1,y2-y1+1);
mag = [y(x1:x2,y1:y2);BB;TGV(x1:x2,y1:y2);BB;...
PCA(x1:x2,y1:y2);BB;BM3d(x1:x2,y1:y2);BB;...
SAIST(x1:x2,y1:y2);BB;WNNM(x1:x2,y1:y2);BB;...
Prop(x1:x2,y1:y2);BB;x(x1:x2,y1:y2)];
mag = imresize(mag,size(d));
ZZZ = [zz,mag,d*7];
ZZZ2 = [zz,mag];

  imwrite(ZZZ,'rec.png')
%  imwrite(ZZZ2,'rec2.png')
figure, imshow(ZZZ)
%%%%%%%%%%%


