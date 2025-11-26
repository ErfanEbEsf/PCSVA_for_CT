function [im_out, PSNR, SSIM]   =  LASSC_Denoising( par )
time0 = clock;
nim   =   par.nim;
ori_im =   par.I;
b     =   par.win;
[h  w ch] =  size(nim);

N = h-b+1;
M = w-b+1;
r = [1:N];
c = [1:M]; 
disp(sprintf('SNR of the noisy image = %.2f \n', csnr(nim, ori_im, 0, 0) ));

im_out = nim;
lamada = par.w;
nsig = par.nSig;

for iter = 1 : par.K        
    im_out  =   im_out + lamada*(nim - im_out);
%     dif     =   im_out-nim;
%     vd      =   nsig^2-(mean(mean(dif.^2)));
%  
%     if (iter==1)
%         par.nSig  = sqrt(abs(vd));            
%     else
%         par.nSig  = sqrt(abs(vd))*lamada;
%     end    
    
  %   if (mod(iter,6)==0) || (iter==1)
        blk_arr = Block_matching( im_out, par);
  %   end
    X = Im2Patch( im_out, par );
    
    Ys = zeros( size(X) );        
    W =  zeros( size(X) );
    L =   size(blk_arr,2);
    for  i =1 : L
        B   =   X(:, blk_arr(:, i));
        mB  =   repmat(mean( B, 2 ), 1, size(B, 2));
        B   =   B-mB;   
        [Ys(:, blk_arr(:,i)), W(:, blk_arr(:,i)), R(i)] =...
          Low_rank_SSC(double(B), par.c1, nsig, mB, par.k0);
    end
    %R_save(iter,:)=R;
    im_out   =  zeros(h,w);
    im_wei   =  zeros(h,w);
    k        =   0;
    for i  = 1:b
        for j  = 1:b
            k    =  k+1;
            im_out(r-1+i,c-1+j)  =  im_out(r-1+i,c-1+j) + reshape( Ys(k,:)', [N M]);
            im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + reshape( W(k,:)', [N M]);
        end
    end
    im_out  =  im_out./(im_wei+eps);
    
    if isfield(par,'I')
        PSNR      =  csnr( im_out, par.I, 0, 0 );
        SSIM      =  cal_ssim( im_out, par.I, 0, 0 );
        RMSE = sqrt(mean((im_out(:) - par.I(:)).^2));
    end
    
    fprintf( 'Iteration %d nSig =%2.2f, SNR =%2.2f, SSIM =%2.4f, RMSE = %.4f\n', iter, par.nSig, PSNR, SSIM,RMSE );
end
if isfield(par,'I')
   PSNR      =  csnr( im_out, par.I, 0, 0 );
   SSIM      =  cal_ssim( im_out, par.I, 0, 0 );
end
%disp(sprintf('Total elapsed time = %f min\n', (etime(clock,time0)/60) ));
return;

%------------------------------------------------------------------
% Re-weighted SV Thresholding
% Sigma = argmin || Y-U*Sigma*V' ||^2 + tau * || Sigma ||_*
%------------------------------------------------------------------
function  [X, W, r]   =  Low_rank_SSC(Y, c1, nsig, m, k)
%% PCA  
data_noisy = (Y');
C = cov(data_noisy);
[V,L,E] = pcacov(C);
W = V(:,1:k);
YY = data_noisy * W;
%YY = wthresh(YY,'h',0.1);
X = YY * W';
X = X';
% r = 1;
% wei = 1;
wei = 1/(sum(E(1:k))*nsig/100); %This is highly in line with BM3D's rule of weight proportionality with reciprocal of total sample variance
W = wei*ones( size(X) );
X = (X+m)*wei;
r = sum(L(1:k));
%% PCA by Laurens - works.
% data_noisy = (Y');
% [YY,Stru] = compute_mapping(data_noisy,'PCA',k);
% %YY = wthresh(YY,'h',0.1);
% W = Stru.M;
% X = YY * W';
% X = X';
% %wei = 1/(sum(E(1:k))*nsig/100); %This is highly in line with BM3D's rule of weight proportionality with reciprocal of total sample variance
% wei = 1;
% W = wei*ones( size(X) );
% X = (X+m)*wei;
% r = 1;
% % r = sum(L(1:k));
%% PCA - SVD
% data_noisy = (Y');
% C = cov(data_noisy);
% [V,L,E] = pcacov(C);
% W = V(:,1:k);
% YY = data_noisy * W;
% 
% [U0,Sigma0,V0]    =   svd(full(YY),'econ');
% Sigma0            =   diag(Sigma0);
% S                 =   max( Sigma0.^2/size(YY, 2) - nsig^2, 0 );
% thr               =   c1*nsig^2./ ( sqrt(S) + eps );
% S                 =   soft(Sigma0, thr);
% %  S = wthresh(Sigma0,'h',thr);
% r                 =   sum( S>0 ); %number of nonzero elemetns in S
%  
% U                 =   U0(:,1:r);
% V                 =   V0(:,1:r);
% YY                =   U*diag(S(1:r))*V';
% 
% X = YY * W';
% X = X';
% wei = 100/(sum(E(1:k))*nsig^2); %This is highly in line with BM3D's rule of weight proportionality with reciprocal of total sample variance
% W = wei*ones( size(X) );
% X = (X+m)*wei;
% r = sum(L(1:k));

return;

