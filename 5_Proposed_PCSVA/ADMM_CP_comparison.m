clear all, close all, clc
%% Iterative Denoising Rec ADMM
% Maybe MP later
clear all; close all; clc
rng(1)
load ('C:\Users\ebrahimesfahani\Desktop\PCSVA\5_Proposed_PCSVA\SAIST_Denoising_Erfan\Data\TestVal_176_P.mat',...
    'X_test_176')
addpath(genpath('SAIST_Denoising_Erfan'))
Op = 0;
FS = 14;
FS2 = 20;
LW = 4.5;

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
N0 = size(X_test_176,3);
N = 2;
X_test_176_rec_ADMM = zeros(m,n,N);

XX = X_test_176(:,:,randi([1,N0],1,N));
%%XX = X_test_176;

for i = 1:N
x = double(XX(:,:,i));
b = Rad(x+Sig_added_noise*randn(m,n));
u = Rad_adj(b);   
u0 = u;
v = zeros(size(u));
w = zeros(size(u));
% SNR_ADMM = zeros(maxit,1);
% SSIM_ADMM = zeros(maxit,1);
% RMSE_ADMM = zeros(maxit,1);

% Main iterations
for j=1:maxit
    tic;
    u = prox_2(v-(1/rho)*w,u0,alpha,rho,gamma);
    v = Image_LASSC_Denoising_Erfan(u+(1/rho)*w,x,Sig);
    w = w + rho*(u-v);
  
    time(i,j) = toc;
    
    SNR_ADMM(i,j)=csnr(u,x,0,0);
    SSIM_ADMM(i,j)=cal_ssim(u,x,0,0);
    RMSE_ADMM(i,j) = sqrt(mean((u(:) - x(:)).^2));
    fprintf('image= %d iteration= %d SNR= %.2f SSIM= %.4f RMSE=%.4f\n',...
    i,j,SNR_ADMM(j),SSIM_ADMM(j),RMSE_ADMM(j));
 end
 tmp_ADMM(:,:,i) = u;
 X_test_176_rec_ADMM(:,:,i) = rgb2gray(insertText(u,[1,1],['ADMM',newline, 'SSIM= ',num2str(100*SSIM_ADMM(i,end),'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_ADMM(i,end),'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));

end

%% Iterative Denoising Rec CP
sigma = 1/10;
tau = 1;
lambda = 1e-5;
n = 176;
m = n;
prox2_sigma = @(r,lambda) r/(lambda*sigma+1);
%% solution by denoising reconstruction

X_test_176_rec_CP = zeros(m,n,N);

for i = 1:N
x = double(XX(:,:,i));
b = Rad(x+Sig_added_noise*randn(m,n));
u = Rad_adj(b);   
u_tild = zeros(m,n);
v_tild = zeros(m,n,2);
r = zeros(size(b));
% SNR_CP = zeros(maxit,1);
% SSIM_CP = zeros(maxit,1);
% RMSE_CP = zeros(maxit,1);

% Main iterations
for j=1:maxit
    r = prox2_sigma (r + sigma*(Rad(u_tild) -b),lambda);
    u_old = u;
    u = u + tau*(- Rad_adj(r));
    tic;
    u = Image_LASSC_Denoising_Erfan(u,x,Sig);
    time(i,j) = toc;
    u_tild = 2*u - u_old;
    
    SNR_CP(i,j)=csnr(u,x,0,0);
    SSIM_CP(i,j)=cal_ssim(u,x,0,0);
    RMSE_CP(i,j) = sqrt(mean((u(:) - x(:)).^2));
    fprintf('image= %d iteration= %d SNR= %.2f SSIM= %.4f RMSE=%.4f\n',...
    i,j,SNR_CP(j),SSIM_CP(j),RMSE_CP(j));
 end
 tmp_CP(:,:,i) = u;
X_test_176_rec_CP(:,:,i) = rgb2gray(insertText(u,[1,1],['CP',newline, 'SSIM= ',num2str(100*SSIM_CP(i,end),'%.2f'),'%',newline,'RMSE= ',num2str(100*RMSE_CP(i,end),'%.2f'),'%'],'BoxOpacity',0,'FontSize',FS,'TextColor','white'));

end

 %% Geet the images side by side
%  for i = 1:N
% z(:,:,i) = [X_test_176_rec_CP(:,:,i),X_test_176_rec_ADMM(:,:,i)];
%  figure, imshow(z(:,:,i))
%  imwrite(z(:,:,i),['ADMM_CP',num2str(i),'.png'])
%  end
 %% Plot convergence curves
for i = 1:N
figHandler= figure;
plot(1:maxit,RMSE_CP(i,:),LineWidth=LW,Color='red')
hold on
plot(1:maxit,RMSE_ADMM(i,:),LineWidth=LW,Color='blue')
legend('CP','ADMM')
set(gca,'FontSize',FS2,'fontWeight','bold')
title('Convergence plots in RMSE')
xlabel('Iteration number')
%ylabel('RMSE')
xticks([1,50,100])
yticks([0, 0.1, 0.2])
save(['ADMM_CP',num2str(i),'.png'],"z")
magnifyOnFigure(...
        figHandler,...
        'units', 'pixels',...
        'magnifierShape', 'rectangle',...
        'initialPositionSecondaryAxes', [250 130 150 150],...
        'initialPositionMagnifier',     [400 92 15 15],...    
        'mode', 'interactive',...    
        'displayLinkStyle', 'straight',...        
        'edgeWidth', 2,...
        'edgeColor', 'black',...
        'secondaryAxesFaceColor', [0.95 0.95 0.95]... 
            ); 
end