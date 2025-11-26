 % Codes by Erfan E. Esfahani.

function [u] = TGV_den(y)
%% Setting parameters
global alpha0 alpha1 lambda b m n D;
maxit = 500;
sigmaTGV = 1/sqrt(12);
tauTGV = sigmaTGV;
lambda = 0.08;
alpha0 = 2;
alpha1 = 1;
[m,n] = size(y);
%% Loadibg data
b = y;

%% Defining operators
prox1 = @(u) (lambda*u + tauTGV*b)/(lambda+tauTGV);
D = @(u) cat(3,dxp(u),dyp(u));
div_1 = @(p) dxm(p(:,:,1)) + dym(p(:,:,2)); 


%% solution by TGV (TGV-MRI)
u = b;        
u_tild = zeros(m,n);
v = D(u);
v_tild = zeros(m,n,2);
p = zeros(m,n,2);
q = zeros(m,n,3);

counter=0;

% Main iterations
%tic;
 for j=1:maxit
    counter = counter + 1;
    p = projP(p + sigmaTGV*(D(u_tild)-v_tild),alpha1);
    q = projQ(q + sigmaTGV*E(v_tild),alpha0);
    u_old = u;
    u = prox1(u + tauTGV*(div_1(p)));
    u_tild = 2*u - u_old;
    v_old = v;
    v = v + tauTGV*(p + div_2(q));
    v_tild = 2*v - v_old;

    fprintf('iteration= %d\n',counter);
 end
%time_TGV = toc;

%% Outputs
%out.IterationsCount = counter;
%out.TGV_Runtime = time_TGV;

function z = E(p)
global m n
z = zeros(m,n,3);
z(:,:,1) = dxm(p(:,:,1));
z(:,:,2) = dym(p(:,:,2));
z(:,:,3) = (dym(p(:,:,1)) + dxm(p(:,:,2)))/2;

function r = div_2(z)
global m n
r = zeros(m,n,2);
r(:,:,1) = dxp(z(:,:,1)) + dyp(z(:,:,3));
r(:,:,2) = dxp(z(:,:,3)) + dyp(z(:,:,2));


function p = projP(p,alpha1)

  absp = sqrt(abs(p(:,:,1)).^2 + abs(p(:,:,2)).^2);
  denom = max(1,absp/alpha1);
  p(:,:,1) = p(:,:,1)./denom;
  p(:,:,2) = p(:,:,2)./denom;  

function q = projQ(q,alpha0)
  absq = sqrt(abs(q(:,:,1)).^2 + abs(q(:,:,2)).^2 + 2*abs(q(:,:,3)).^2);
  denom = max(1,absq/alpha0);
  q(:,:,1) = q(:,:,1)./denom;
  q(:,:,2) = q(:,:,2)./denom;
  q(:,:,3) = q(:,:,3)./denom;  
 
             