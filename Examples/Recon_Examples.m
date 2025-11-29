function Recon_Examples
clear all, close all, clc
load Patient16.mat
load PAT016.mat
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\8_Examples')
i = 70; %148
x = X(:,:,i); mm = Mask(:,:,i);
x = imresize(x,176/512);
mm = imresize(mm,176/512);
% mm = bwperim(mm,4);
% mm = boundarymask(mm);

% % bw = activecontour(x*255,mm,1000,'edge');
% % imshow(x)
% % hold on
% % h = visboundaries(bw,'Color','b')
% % hold off
% % h.Children(1)
% % h.Children(2)
% % h.Children(2).Visible = 'off';
% % h.Children(1).LineWidth = 1;
% 
 x3 = repmat(x,1,1,3);
 Seg = double(labeloverlay(x,mm,'Transparency',0 ));
 Seg = Seg/max(Seg(:));
% %figure, imshow([x3, Seg])
%  figure, imshow([x3,Seg])
%% Setting parameters 
maxit = 100;
maxitTGV = 500; %500
sigma = 1/10;
tau = 1;
lambda = 1e-5;
n = 176; % size(x,1);
m = n;
Sig = 0.03; %.04
Sig_added_noise = 0;
x = imresize(x,[m,n]);
%% Defining operators
theta = 0:5*1:179; % *2
Rad =@(x) radon(x,theta);
Rad_adj = @(r) iradon(r,theta(2)-theta(1),n);
prox2_sigma = @(r,lambda) r/(lambda*sigma+1);

%% Proposed Wiener
% cd('C:\Users\ebrahimesfahani\Desktop\5_Proposed_PCA_Wiener')
% addpath(genpath('Utils'))
% tic,
% b = Rad(x+Sig_added_noise*randn(m,n));
% u = Rad_adj(b);   
% u_tild = zeros(m,n);
% v_tild = zeros(m,n,2);
% r = zeros(size(b));
% SNR_Proposed = zeros(maxit,1);
% SSIM_Proposed = zeros(maxit,1);
% RMSE_Proposed = zeros(maxit,1);
% 
% % Main iterations
% for j=1:maxit
%     r = prox2_sigma (r + sigma*(Rad(u_tild) -b),lambda);
%     u_old = u;
%     u = u + tau*(- Rad_adj(r));   
%     u = Image_LASSC_Denoising_Erfan(u,x,Sig);
%     u_tild = 2*u - u_old;
%     
%     SNR_Proposed(j)=csnr(u,x,0,0);
%     SSIM_Proposed(j)=cal_ssim(u,x,0,0);
%     RMSE_Proposed(j) = sqrt(mean((u(:) - x(:)).^2));
%     fprintf('image= %d iteration= %d SNR= %.2f SSIM= %.4f RMSE=%.4f\n',...
%     i,j,SNR_Proposed(j),SSIM_Proposed(j),RMSE_Proposed(j));
% end
% Prop = u;
% SNR_Prop = SNR_Proposed(end);SSIM_Prop = SSIM_Proposed(end);
% RMSE_Prop = RMSE_Proposed(end);
% time = toc

%% Proposed
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\5_Proposed_PCSVA')
addpath(genpath('SAIST_Denoising_Erfan'))
tic,
b = Rad(x+Sig_added_noise*randn(m,n));
u = Rad_adj(b);  
fbp = u;
u_tild = zeros(m,n);
r = zeros(size(b));
SNR_Proposed = zeros(maxit,1);
SSIM_Proposed = zeros(maxit,1);
RMSE_Proposed = zeros(maxit,1);

% Main iterations
for j=1:maxit
    r = prox2_sigma (r + sigma*(Rad(u_tild) -b),lambda);
    u_old = u;
    u = u + tau*(- Rad_adj(r));   
    u = Image_LASSC_Denoising_Erfan(u,x,Sig);
    u_tild = 2*u - u_old;
    
    SNR_Proposed(j)=csnr(u,x,0,0);
    SSIM_Proposed(j)=cal_ssim(u,x,0,0);
    RMSE_Proposed(j) = sqrt(mean((u(:) - x(:)).^2));
    fprintf('image= %d iteration= %d SNR= %.2f SSIM= %.4f RMSE=%.4f\n',...
    i,j,SNR_Proposed(j),SSIM_Proposed(j),RMSE_Proposed(j));
end
Prop = u;
SNR_Prop = SNR_Proposed(end);SSIM_Prop = SSIM_Proposed(end);
RMSE_Prop = RMSE_Proposed(end);
time = toc

%% WNNM
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\4_WNNM_2')
addpath(genpath('SAIST_Denoising_Erfan'))
tic,
b = Rad(x+Sig_added_noise*randn(m,n));
u = Rad_adj(b);   
u_tild = zeros(m,n);
r = zeros(size(b));
SNR_Proposed = zeros(maxit,1);
SSIM_Proposed = zeros(maxit,1);
RMSE_Proposed = zeros(maxit,1);

% Main iterations
for j=1:maxit
    r = prox2_sigma (r + sigma*(Rad(u_tild) -b),lambda);
    u_old = u;
    u = u + tau*(- Rad_adj(r));   
    u = Image_LASSC_Denoising_Erfan(u,x,Sig);
    u_tild = 2*u - u_old;
    
    SNR_Proposed(j)=csnr(u,x,0,0);
    SSIM_Proposed(j)=cal_ssim(u,x,0,0);
    RMSE_Proposed(j) = sqrt(mean((u(:) - x(:)).^2));
    fprintf('image= %d iteration= %d SNR= %.2f SSIM= %.4f RMSE=%.4f\n',...
    i,j,SNR_Proposed(j),SSIM_Proposed(j),RMSE_Proposed(j));
end
WNNM = u;
SNR_WNNM = SNR_Proposed(end);SSIM_WNNM = SSIM_Proposed(end);
RMSE_WNNM = RMSE_Proposed(end);
time = toc

%% SAIST
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\3_SAIST')
addpath(genpath('SAIST_Denoising_Erfan'))
tic,
b = Rad(x+Sig_added_noise*randn(m,n));
u = Rad_adj(b);   
u_tild = zeros(m,n);
r = zeros(size(b));
SNR_Proposed = zeros(maxit,1);
SSIM_Proposed = zeros(maxit,1);
RMSE_Proposed = zeros(maxit,1);

% Main iterations
for j=1:maxit
    r = prox2_sigma (r + sigma*(Rad(u_tild) -b),lambda);
    u_old = u;
    u = u + tau*(- Rad_adj(r));   
    u = Image_LASSC_Denoising_Erfan(u,x,Sig);
    u_tild = 2*u - u_old;
    
    SNR_Proposed(j)=csnr(u,x,0,0);
    SSIM_Proposed(j)=cal_ssim(u,x,0,0);
    RMSE_Proposed(j) = sqrt(mean((u(:) - x(:)).^2));
    fprintf('image= %d iteration= %d SNR= %.2f SSIM= %.4f RMSE=%.4f\n',...
    i,j,SNR_Proposed(j),SSIM_Proposed(j),RMSE_Proposed(j));
end
SAIST = u;
SNR_SAIST = SNR_Proposed(end);SSIM_SAIST = SSIM_Proposed(end);
RMSE_SAIST = RMSE_Proposed(end);
time = toc

%% BM3D
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\2_BM3D')
tic,
b = Rad(x+Sig_added_noise*randn(m,n));
u = Rad_adj(b);   
u_tild = zeros(m,n);
r = zeros(size(b));
SNR_Proposed = zeros(maxit,1);
SSIM_Proposed = zeros(maxit,1);
RMSE_Proposed = zeros(maxit,1);

% Main iterations
for j=1:maxit
    r = prox2_sigma (r + sigma*(Rad(u_tild) -b),lambda);
    u_old = u;
    u = u + tau*(- Rad_adj(r));   
    [~,u] = BM3D2(1,u,255*Sig);
    u = double(u);
    u_tild = 2*u - u_old;
    
    SNR_Proposed(j)=csnr(u,x,0,0);
    SSIM_Proposed(j)=cal_ssim(u,x,0,0);
    RMSE_Proposed(j) = sqrt(mean((u(:) - x(:)).^2));
    fprintf('image= %d iteration= %d SNR= %.2f SSIM= %.4f RMSE=%.4f\n',...
    i,j,SNR_Proposed(j),SSIM_Proposed(j),RMSE_Proposed(j));
end
BM3d = u;
SNR_BM3d = SNR_Proposed(end);SSIM_BM3d = SSIM_Proposed(end);
RMSE_BM3d = RMSE_Proposed(end);
time = toc

%% TGV
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\1_TGV')
addpath(genpath('./UTILS'));
sigmaTGV = 1/sqrt(12);
tauTGV = 1/sqrt(12);
alpha0 = 2;
alpha1 = 0.01;

D = @(u) cat(3,dxp(u),dyp(u));
div_1 = @(p) dxm(p(:,:,1)) + dym(p(:,:,2)); 

tic,
b = Rad(x+Sig_added_noise*randn(m,n));
u = Rad_adj(b);   
u_tild = zeros(m,n);
v = D(u);
v_tild = zeros(m,n,2);
p = zeros(m,n,2);
q = zeros(m,n,3);
r = zeros(size(b));
counter=0;

SNR_Proposed = zeros(maxit,1);
SSIM_Proposed = zeros(maxit,1);
RMSE_Proposed = zeros(maxit,1);

% Main iterations
for j=1:maxitTGV
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
TGV = u;
SNR_TGV = SNR_TGV(end);SSIM_TGV = SSIM_TGV(end);
RMSE_TGV = RMSE_TGV(end);
time = toc

%% Classical PCA
cd('C:\Users\ebrahimesfahani\Desktop\PCSVA\6_classical_PCA')
addpath(genpath('SAIST_Denoising_Erfan'))
tic,
b = Rad(x+Sig_added_noise*randn(m,n));
u = Rad_adj(b);   
u_tild = zeros(m,n);
r = zeros(size(b));
SNR_Proposed = zeros(maxit,1);
SSIM_Proposed = zeros(maxit,1);
RMSE_Proposed = zeros(maxit,1);

% Main iterations
for j=1:maxit
    r = prox2_sigma (r + sigma*(Rad(u_tild) -b),lambda);
    u_old = u;
    u = u + tau*(- Rad_adj(r));   
    u = Image_LASSC_Denoising_Erfan(u,x,Sig);
    u_tild = 2*u - u_old;
    
    SNR_Proposed(j)=csnr(u,x,0,0);
    SSIM_Proposed(j)=cal_ssim(u,x,0,0);
    RMSE_Proposed(j) = sqrt(mean((u(:) - x(:)).^2));
    fprintf('image= %d iteration= %d SNR= %.2f SSIM= %.4f RMSE=%.4f\n',...
    i,j,SNR_Proposed(j),SSIM_Proposed(j),RMSE_Proposed(j));
end
PCA = u;
SNR_PCA = SNR_Proposed(end);SSIM_PCA = SSIM_Proposed(end);
RMSE_PCA = RMSE_Proposed(end);
time = toc
%% Output
SSIM_FBP = cal_ssim(fbp,x,0,0);
RMSE_FBP = sqrt(mean((fbp(:) - x(:)).^2));
Delta = 50;
x1 = 100; x2 = x1 + Delta; %70 %120
y1 = 120; y2 = y1 + Delta; %78 %117 
Op = 0;
FS = 16;

BB = zeros(1,size(fbp,2));
z = [fbp;BB;TGV;BB;PCA;BB;BM3d;BB;...
    SAIST;BB;WNNM;BB;Prop;BB;x];
d = abs(z - [x;BB;x;BB;x;BB;x;BB;...
    x;BB;x;BB;x;BB;x]);

x = rgb2gray(insertText(x,[1,1],['Ground Truth',newline,'SSIM= 100%',newline, 'RMSE= 0%'],'BoxOpacity',Op,'FontSize',FS,TextColor='white'));
fbp = rgb2gray(insertText(fbp,[1,1],['FBP',newline, 'SSIM= ',num2str(100*SSIM_FBP,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_FBP,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));
PCA = rgb2gray(insertText(PCA,[1,1],['PCA',newline, 'SSIM= ',num2str(100*SSIM_PCA,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_PCA,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));
TGV = rgb2gray(insertText(TGV,[1,1],['TGV',newline, 'SSIM= ',num2str(100*SSIM_TGV,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_TGV,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white')); 
SAIST = rgb2gray(insertText(SAIST,[1,1],['SAIST',newline, 'SSIM= ',num2str(100*SSIM_SAIST,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_SAIST,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));
BM3d = rgb2gray(insertText(BM3d,[1,1],['BM3D',newline, 'SSIM= ',num2str(100*SSIM_BM3d,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_BM3d,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));
WNNM = rgb2gray(insertText(WNNM,[1,1],['WNNM',newline, 'SSIM= ',num2str(100*SSIM_WNNM,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_WNNM,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));
Prop = rgb2gray(insertText(Prop,[1,1],['Proposed',newline, 'SSIM= ',num2str(100*SSIM_Prop,'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_Prop,'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));

zz = [fbp;BB;TGV;BB;PCA;BB;BM3d;BB;...
    SAIST;BB;WNNM;BB;Prop;BB;x];


BB = zeros(1,y2-y1+1);
mag = [fbp(x1:x2,y1:y2);BB;TGV(x1:x2,y1:y2);BB;...
PCA(x1:x2,y1:y2);BB;BM3d(x1:x2,y1:y2);BB;...
SAIST(x1:x2,y1:y2);BB;WNNM(x1:x2,y1:y2);BB;...
Prop(x1:x2,y1:y2);BB;x(x1:x2,y1:y2)];
mag = imresize(mag,size(d));
ZZZ = [zz,mag,d*10];
ZZZ2 = [zz,mag];

%imwrite(ZZZ,'C:\Users\ebrahimesfahani\Desktop\PCSVA\8_Examples\rec.png')
%imwrite(ZZZ2,'C:\Users\ebrahimesfahani\Desktop\PCSVA\8_Examples\rec2.png')
figure, imshow(ZZZ)

%%

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