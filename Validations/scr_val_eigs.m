
% This script (1) checks that the eigenvalues are the same for the three
% separate uncoupled governing equations for mathcalD, mathcalR, tilp; and (2)
% that these eigenvalues agree with those in Ioannou and Lindzen (1993) (the
% forcing and parameters are chosen to match those of the the latter study).
%
% Notes:
%
% 1) Below, we make use of the fact that the governing matrix operators (Ltilp,
% LtilmfD, LtilmfR) show uncoupling between the even and odd harmonics (note the
% sparsity pattern whereby a row either couples terms of degree even or odd but
% not both.) This allows a separation of the matrices into two (symmetric and an
% asymmetric) as below.
% 
% 2) The example here uses "symmetric" nominally to refer to modes with the same
% degree parity as s. When s is even (e.g. in s=2 example below) "symmetric"
% means the eigenfunctions are symmetric with respect to equator.   
%

%% Housekeeping:

clear
TROPFPATH  = '~/Desktop/TROPF';     % set path of TROPF root directory
cd(TROPFPATH); scr_startup_tropf;   % run startup script to set paths to TROPF libraries
workdir    = [TROPFPATH,'/Validations'];
cd(workdir)



%% % Select example input parameters:

% We consider different cases; the parameters that differ are collected here up front (so just these change) between cases:

N     = 500 ; % number of terms in SH expansion


tilalpd = 0; tilalpr = 0; tilalpp = 0; % inviscid following IL1993 assumptions 
%tilalpd = rand(1); tilalpr = rand(1); tilalpp = rand(1); % test TROPF equations which allow for dissipation and complex eigenvalues



%%

% initialize:
% method parameters:
tilOm = 1    ; % nondim rotation rate (\Omega / \Omega_s)
%N     = 1500 ; % number of terms in SH expansion


% forcing parameters:
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1); % intitialize

 
% Case sectorial force (P22) with Io-on-Jupiter tilom:
nF = 2; sF = 2; s = sF; % degree and order of forcing potential 
Ntrunc = N + sF - 1;    % truncation degree
PnFsF_amp  = 3;         % amplitude of associated Legendre function of degree nF, order sF;
tilom      = -0.766  ;  % frequency of forcing potential 
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = 1i*( (1/2) / PnFsF_amp ); % tilmfG normalized to have unit amplitude
tilcesq = 1; % this choice winds up being just a reference value 
%tilalpd = 0; tilalpr = 0; tilalpp = 0; % inviscid following IL1993 assumptions 
% prescribe squared slowness coeff vector:
tilnusqns =  (1 + 1i*tilalpp/tilom)./tilcesq;
% dissipation and slowness operators:
nvec      =  nVec(N,s);
Lalphad   =  build_Lalphad(nvec, tilalpd);
Lalphar   =  build_Lalphar(nvec, tilalpr);
LV        =  build_LV_fromSlowness(N,tilnusqns);
 

% Other operators:
LL        =  build_LL(nvec);
LA        =  build_LA(nvec,tilOm, s,tilom, Lalphad, LL);
LB        =  build_LB(nvec,tilOm, s,tilom, Lalphar, LL);
LC        =  build_LC(nvec,tilOm, s);
LD        =  build_LD(nvec,tilOm, s,tilom, Lalphad,LV, LL);
LVi       =  build_LVi_fromSlowness(N,tilnusqns);
LLi       =  build_LLi(nvec);
LAi       =  build_LAi(nvec,tilOm, s,tilom, Lalphad, LL);
LBi       =  build_LBi(nvec,tilOm, s,tilom, Lalphar,LL);
LDi       =  build_LDi(nvec,tilOm, s,tilom, Lalphad,LV, LL);

Ltilp     = build_Ltilp(tilom, LV, LLi,LA,LC,LBi);
LtilmfD   = build_LtilmfD(LLi,LC,LBi,LD);
LtilmfR   = build_LtilmfR(N,tilOm,  Lalphad,Lalphar,LV,  tilom,s);



%% Eigenvalues in Ltilp formulation:
 
LVjunk = - ( LLi*(LA - LC*LBi*LC)*tilom*LLi ) \ spdiags(ones(N,1), 0,N,N);
LVjunk_sym = LVjunk(1:2:end,1:2:end); % matrix for symmetric 
LVjunk_asy = LVjunk(2:2:end,2:2:end); % matrix for asymmetric
[V,D]  = eig(full(LVjunk_sym)); 
tilcesqvals  = 1./diag(D); % squared eigen-wavespeeds
[junk,ii_descend] = sort(abs(tilcesqvals),'descend');
tilcesqvals_sym_tilp    = tilcesqvals(ii_descend,:);
[V,D]  = eig(full(LVjunk_asy)); 
tilcesqvals  = 1./diag(D); % squared eigen-wavespeeds
[junk,ii_descend] = sort(abs(tilcesqvals),'descend');
tilcesqvals_asy_tilp    = tilcesqvals(ii_descend,:);


%% Eigenvalues in LtilmfD formulation:
 
LVjunk = - tilom*( LLi*(LA - LC*LBi*LC)*LLi ) 
LVjunk_sym = LVjunk(1:2:end,1:2:end); % matrix for symmetric 
LVjunk_asy = LVjunk(2:2:end,2:2:end); % matrix for asymmetric
[V,D]  = eig(full(LVjunk_sym)); 
tilcesqvals  = diag(D); % squared eigen-wavespeeds
[junk,ii_descend] = sort(abs(tilcesqvals),'descend');
tilcesqvals_sym_tilmfD    = tilcesqvals(ii_descend,:);
[V,D]  = eig(full(LVjunk_asy)); 
tilcesqvals  = diag(D); % squared eigen-wavespeeds
[junk,ii_descend] = sort(abs(tilcesqvals),'descend');
tilcesqvals_asy_tilmfD    = tilcesqvals(ii_descend,:);




%% Eigenvalues in LtilmfR formulation:
 
LVjunk =  tilom*LLi* (LC*LBi -LA*inv(LC)) * LC*LLi;
LVjunk_sym = LVjunk(1:2:end,1:2:end); % matrix for symmetric 
LVjunk_asy = LVjunk(2:2:end,2:2:end); % matrix for asymmetric
[V,D]  = eig(full(LVjunk_sym)); 
tilcesqvals  = diag(D); % squared eigen-wavespeeds
[junk,ii_descend] = sort(abs(tilcesqvals),'descend');
tilcesqvals_sym_tilmfR    = tilcesqvals(ii_descend,:);
[V,D]  = eig(full(LVjunk_asy)); 
tilcesqvals  = diag(D); % squared eigen-wavespeeds
[junk,ii_descend] = sort(abs(tilcesqvals),'descend');
tilcesqvals_asy_tilmfR    = tilcesqvals(ii_descend,:);



%% For consistency check, show the eigenvalues are the same for the three different formulations:


enum = 1:N;
 
figure(11)
plot(...
     enum(1:2:end),tilcesqvals_sym_tilp,'rhex',    enum(2:2:end),tilcesqvals_asy_tilp,'bsquare',  ...
     enum(1:2:end),tilcesqvals_sym_tilmfD,'rdiamond',  enum(2:2:end),tilcesqvals_asy_tilmfD,'bpent',  ...
     enum(1:2:end),tilcesqvals_sym_tilmfR,'r+',  enum(2:2:end),tilcesqvals_asy_tilmfR,'b.',  ...
    'MarkerSize',15) 
if mod(s,2), legend('A asymmetric','A symmetric','B asymmetric','B symmetric','C asymmetric','C symmetric'), else, legend('A symmetric','A asymmetric','B symmetric','B asymmetric','C symmetric','C asymmetric'), end
set(gca,'xlim',[0 30])
ylabel('eigenvalues (i.e. ${\tilde c_{e,(k)}^2}$)', 'interpreter','latex')
grid

%print -dpng eigenvaluesIL
print -dpng eigenvaluesWTropfDissip


if imag(tilcesqvals),
     
figure(12)
plot(...
     enum(1:2:end),imag(tilcesqvals_sym_tilp),'rhex',    enum(2:2:end),imag(tilcesqvals_asy_tilp),'bsquare',  ...
     enum(1:2:end),imag(tilcesqvals_sym_tilmfD),'rdiamond',  enum(2:2:end),imag(tilcesqvals_asy_tilmfD),'bpent',  ...
     enum(1:2:end),imag(tilcesqvals_sym_tilmfR),'r+',  enum(2:2:end),imag(tilcesqvals_asy_tilmfR),'b.',  ...
    'MarkerSize',15) 
if mod(s,2), legend('A asymmetric','A symmetric','B asymmetric','B symmetric','C asymmetric','C symmetric'), else, legend('A symmetric','A asymmetric','B symmetric','B asymmetric','C symmetric','C asymmetric'), end
set(gca,'xlim',[0 30])
ylabel('eigenvalues (i.e. ${\tilde c_{e,(k)}^2}$)', 'interpreter','latex')
grid
print -dpng eigenvaluesWTropfDissipImag
end



%% Make Table


junk = reshape(tilcesqvals_sym_tilp(1:60),20,3)
dlmwrite('symN120',junk,'delimiter','\t','precision','%10.7e')


junk = reshape(tilcesqvals_asy_tilp(1:60),20,3)
dlmwrite('asyN120',junk,'delimiter','\t','precision','%10.7e')




   





