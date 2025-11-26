function variance_estimate=iedd(im)

% The function estimates noise variance in the grayscale image im
%
% Please refer to: 
% M.Ponomarenko, N.Gapon, V.Voronin, K. Egiazarian, "Blind estimation 
% of white Gaussian noise variance in highly textured images",
% Image Processing: Algorithms and Systems XVI, 2018, 5p.

  blks=dctm(im);
  ene = blks.^2; ene = sum(ene');
  [m1,m2] = sort(ene);
  K=49;
  pz=find(m2==64);
  m2=m2(1:K);
  if (pz<=K)&&(m1(pz)<m1(1)*1.3)
    m2(pz)=m2(1);m2(1)=64;
  end
  m = mymad(blks(m2(1),:));
  for it=1:3
     z = blks(m2,:);
     y=mean(z.^2);
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
  T=dct(eye(blockSize));
  blks = kron(T,T)*blks;
  [y x]=size(blks);
end

function r=mymad(d)
  d=d(:);
  m=median(d);
  r=median(abs(d-m))*1.4826;
end