%%
clear all; close all; clc
rng('default')
load ('C:\Users\ebrahimesfahani\Desktop\PCSVA\5_Proposed_PCSVA\SAIST_Denoising_Erfan\Data\TestVal_176_P.mat',...
    'X_test_176')
addpath(genpath('SAIST_Denoising_Erfan')) 
%%
N = length(X_test_176);

parfor i = 1:N
 GaussSig(i) = sqrt(iedd(X_test_176(:,:,i)));
end
min_GaussSig = min(GaussSig)
max_GaussSig = max(GaussSig)

%% 