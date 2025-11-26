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
          Low_rank_SSC( double(B), par.C, nsig,mB);
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
% [~,im_out] = BM3D_Wiener(im_out,nim,nsig*255);
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
function  [X, W, r]   =   Low_rank_SSC( Y, C, nsig, m )
%% WNNM 1 - failed
% [U0,Sigma0,V0]    =   svd(full(Y),'econ');
% Sigma0            =   diag(Sigma0);
% S                 =   max( Sigma0.^2- size(Y,2)*nsig^2, 0);
% thr               =   (c1*sqrt(size(Y,2)))./ ( sqrt(S) + eps );
% %mean(thr)
% S                 =   soft(Sigma0, thr);
% r                 =   sum( S>0 ); %number of nonzero elemetns in S
% 
% U                 =   U0(:,1:r);
% V                 =   V0(:,1:r);
% X                 =   U*diag(S(1:r))*V';
% wei = 1; 
% % if r==size(Y,1) % Weights not given in WNNM paper
% %     wei           =   1/size(Y,1);
% % else
% %     wei           =   (size(Y,1)-r)/size(Y,1);
% %end
% W                 =   wei*ones( size(X) );
% X                 =   (X+m)*wei;

 %% WNNM 2
% [U0,Sigma0,V0]    =   svd(full(Y),'econ');
% Sigma0            =   diag(Sigma0);
% S = zeros(size(Sigma0));
% C1 = Sigma0 - eps; C2 = (Sigma0 + eps).^2 - 4*C;
% for i = 1:length(Sigma0)
%     if C2(i)<0
%         S(i) = 0;
%     else 
%         S(i) = 0.5*(C1(i) + sqrt(C2(i)));
%     end
% end
% 
% r                 =   sum( S>0 );
% U                 =   U0;
% V                 =   V0;
% X                 =   U*diag(S)*V';
% 
% if r==size(Y,1) % Weights not given in WNNM paper
%     wei           =   1/size(Y,1);
% else
%     wei           =   (size(Y,1)-r)/size(Y,1);
% end
% W                 =   wei*ones( size(X) );
% X                 =   (X+m)*wei;
%% WNNM 3
 [U,SigmaY,V] =   svd(full(Y),'econ');    
    PatNum       = size(Y,2);
    TempC  = C*sqrt(PatNum)*2*nsig^2;
    temp=(SigmaY-eps).^2-4*(TempC-eps*SigmaY);
    ind=find (temp>0);
    svp=length(ind);
    SigmaX=max(SigmaY(ind)-eps+sqrt(temp(ind)),0)/2;
    X =  U(:,1:svp)*diag(SigmaX)*V(:,1:svp)' + m; 
    wei = 1;
    W = wei*ones( size(X) );
    r = svp;
return;

