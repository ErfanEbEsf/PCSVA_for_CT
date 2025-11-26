clear all
close all
clc
rng('default')
addpath(genpath('Utils'))
%% 
n = 176;
N_MC = 100;
U = randn([n,n,N_MC]);
U = max(0,U);
for k = 1:N_MC
U(:,:,k) = U(:,:,k)./max(max(U(:,:,k)));
end
%%
Sig = [0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01];
j = 1;
E = zeros(N_MC,length(Sig));
for sig = Sig
parfor i = 1:N_MC
U_den = Image_LASSC_Denoising_Erfan...
 (U(:,:,i),U(:,:,i),sig);
E(i,j) = norm(U_den-U(:,:,i));
end
j = j+1;
end
%% 
figure,
plot(Sig,max(E,[],1),LineWidth=3,Marker="*")
xlabel('\sigma')
ylabel('max_u ||H_\sigma(u) - u||_2')
grid on
set(gca,'FontSize',12) 