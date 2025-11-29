function variance_estimate=iedd2(im)
  blks=dctm(im);
  ene = blks.^2; ene = sum(ene'); %energy of each DCT coeff
  [m1,m2] = sort(ene);  % sort by ascending energy
  K=49;
  pz=find(m2==64);%Where the last DCT coeff ranks in the sorted vec
  m2=m2(1:K);  %K lowst-energy DCT coeffs are seleced
%   if (pz<=K)&&(m1(pz)<m1(1)*1.3)
%   % if the last DCT coeff is among the lowest energy coeffs && 
%   % its energy is less than 1.3 times min energy 
%      m2(pz)=m2(1);m2(1)=64;
%     % put the lowest energy coeff where the highers freq coeff ranks
%     % && set the highers freq coeff as the lowest energy one
%   end
  m = mymad(blks(m2(1),:)); %Grab the lowst energy DCT coeff and do the median thing on it
  for it=1:1
     z = blks(m2,:); %K lowst-energy DCT coeffs are seleced
     y=mean(z.^2);%(local?) varriance for each of the lowest energy coeffs are calculated
     map=y<(1+sqrt(8/K))*m^2;
     if nnz(map)>1024
       z=z(:,map); m=mymad(z(1,:));
     end
  end
  variance_estimate = m^2;
end

function blks = dctm(im)
  blockSize=8;
  blks = reshape(im2col(double(im), [blockSize blockSize],'sliding'), blockSize*blockSize, []);
  T=mirt_dctn(eye(blockSize));
  blks = kron(T,T)*blks;
  [y x]=size(blks);
end

function r=mymad(d)
  d=d(:);
  m=median(d);
  r=median(abs(d-m))*1.4826;
end