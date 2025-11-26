clear all, close all
x1 = imread("fig1_phantom.png");
x2 = imread("fig1_phantom_zoom.png");
x3 = imread("fig1_patch_noisy.png");
x4 = imread("fig1_patch_PCA.png");
x5 = imread("fig1_patch_prop.png");
x6 = imread("fig1_patch_OG.png");
N = 512;
N_v = 100;

x1 = imresize(x1,[N,N]);
x2 = imresize(x2,[N,N]);
x3 = imresize(x3,[N,N]);
x3 = repmat(x3,1,1,3);
x4 = imresize(x4,[N,N]);
x4 = repmat(x4,1,1,3);
x5 = imresize(x5,[N,N]);
x5 = repmat(x5,1,1,3);
x6 = imresize(x6,[N,N]);
x6 = repmat(x6,1,1,3);
%%
z_h = 255*ones(size(x1,1),20,3);
X1 = cat(2,x1,z_h,x2,z_h,x3);
z_v = 255*ones(N_v,size(X1,2),3);
X2 = cat(2,x4,z_h,x5,z_h,x6);
X = cat(1,X1,z_v,X2);
%%
% X1 = cat(2,x1,z_h,x2);
% z_v = 255*ones(150,size(X1,2),3);
% 
% X2 = cat(2,x3,z_h,x4);
% X3 = cat(2,x5,z_h,x6);
% X = cat(1,X1,z_v,X2,z_v,X3,z_v);
figure, imshow(X)