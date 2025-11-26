clear all, close all, clc
%% RMSE
RMSE_sig_5_PCA = [2.59,2.11,1.87,1.79,1.8,...
    1.87,2,2.15];
k_sig_5_PCA = 1:8;

RMSE_sig_10_PCA = [3.09,2.65,2.5,2.5,2.6,...
    2.77,2.99,3.24];
k_sig_10_PCA=1:8;

RMSE_sig_5_Prop =[1.55,1.51,1.5,1.49,1.49,...
    1.5,1.51,1.52,1.53,1.55];
k_sig_5_Prop = 6:2:24;

RMSE_sig_10_Prop=[2.64,2.34,2.29,2.27,2.26,...
    2.27,2.29,2.30,2.30,2.30,2.3];
k_sig_10_Prop=[2:2:4,5,6:2:20];

figure,
plot(k_sig_5_PCA,RMSE_sig_5_PCA,LineWidth=1.5)
hold on
plot(k_sig_5_Prop,RMSE_sig_5_Prop,LineWidth=1.5)
hold on
title('Tuning of k at \sigma_{noise}= 0.05')
xlabel('Number of retained PCs (k)')
ylabel('Average RMSE on tuning set (%)')
legend('PCA','Proposed')
grid on

figure,
plot(k_sig_10_PCA,RMSE_sig_10_PCA,LineWidth=1.5)
hold on
plot(k_sig_10_Prop,RMSE_sig_10_Prop,LineWidth=1.5)
title('Tuning of k at \sigma_{noise}= 0.10')
legend('PCA','Proposed')
xlabel('Number of retained PCs (k)')
ylabel('Average RMSE on tuning set (%)')
grid on
%% SSIM

SSIM_sig_10_PCA=[86.48,86.71,85.11,...
    82.38,78.92];
k_sig_10_PCA_SSIM=[1:1:5];

SSIM_sig_10_Prop=[88.53,89.40,89.47,89.43,89.23...
    88.95,88.69,88.60];
k_sig_10_Prop_SSIM=[2:2:4,5,6:2:14];


figure,
plot(k_sig_10_PCA_SSIM,SSIM_sig_10_PCA,LineWidth=1.5)
hold on
plot(k_sig_10_Prop_SSIM,SSIM_sig_10_Prop,LineWidth=1.5)



title('Tuning of k at \sigma_{noise}= 0.10')
legend('PCA','Proposed')
xlabel('Number of retained PCs (k)')
ylabel('Average SSIM on tuning set (%)')
