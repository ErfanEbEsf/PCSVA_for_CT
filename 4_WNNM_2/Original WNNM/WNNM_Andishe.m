function  [X] =  WNNM_Andishe( Y, C, NSig, m )
    [U,SigmaY,V] =   svd(full(Y),'econ');    
    PatNum       = size(Y,2);
    TempC  = C*sqrt(PatNum)*2*NSig^2;
    temp=(SigmaY-eps).^2-4*(TempC-eps*SigmaY);
    ind=find (temp>0);
    svp=length(ind);
    SigmaX=max(SigmaY(ind)-eps+sqrt(temp(ind)),0)/2;
    X =  U(:,1:svp)*diag(SigmaX)*V(:,1:svp)' + m;     
return;
% 
% [U,Sigma0,V] =   svd(full(Y),'econ');    
%     Sigma0 = diag(Sigma0);
%     PatNum       = size(Y,2);
%     TempC  = C*sqrt(PatNum)*2*NSig^2;
%     S =   max( Sigma0.^2/size(Y, 2) - NSig^2, 0 );
%     thr =   TempC./ ( sqrt(S) + eps );
%     S   =   soft(Sigma0, thr);
%     %temp=(SigmaY-eps).^2-4*(TempC-eps*SigmaY);
%     ind=find (S>0);
%     svp=length(ind);
%     %SigmaX=max(Sigma0(ind)-eps+sqrt(temp(ind)),0)/2;
%     X =  U(:,:)*diag(S)*V(:,:)' + m;     
