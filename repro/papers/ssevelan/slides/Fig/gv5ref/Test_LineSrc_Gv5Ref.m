%% Model
clear all;
startup;
%matlabpool open local 12;
% grid
x = 0:10:3000;
z = 0:10:1500;

[o,d,n] = grid2odn(z,x);

% velocity
[xx,zz]=meshgrid(x,z);

vel = 2000*ones(n);
[xx,zz]=meshgrid(x,z);
w1 = zeros(n(1),1);   i0=find(abs(z-600)<=400);     w1(i0) = tukeywin(length(i0),0.4);
w2 = zeros(n(2),1);   j0=find(abs(x-1500)<=600);    w2(j0) = tukeywin(length(j0),0.2);

gv =  -600*exp(-((xx-1500).^2/300^2+(zz-600).^2/300^2));


vref = zeros(n);
vref(21,:) = 500;
vref(41,:) = 500;
vref(61,:) = 500;
vref(81,:) = 500;
vref(101,:) = 500;

S = opKron(opSmooth(n(2),2),opSmooth(n(1),2));

v = vel + reshape(vref(:),n) + (w1*w2').*gv;
v0 = vel;

%% Smooth operator

Bx = opSpline1D([100:100:2900],x,[0 0]); %% dirichlet bc
Bz = opSpline1D([100:100:1600],z,[0 0]); %% dirichlet bc;
B  = opKron(Bx, Bz);

beta = 0;
mul  = 1;
z1=z/1000;x1=x/1000;
Lap = opFunction(prod(n),prod(n),{@(g)lapsmoother2(g, z1, x1, beta,mul),@(g)lapsmoother2(g, z1, x1, beta,mul)},0,1);


% parameters

% model grid
model.o = o;
model.d = d;
model.n = n;

% absorbing boundary
model.nb = [40 40]; %% the number of sampling in the ABCs domain.

% source/receiver grid
model.xsrc = [500:100:2500];
model.xsrc = 1500;
model.zsrc = 40;
model.xrec = 20:20:2980;
model.zrec = 20;

%% time sample and frequency sampling rate via Nyquist law.
% time 
t  = 0:.004:4;   %% t=0:dt:T;
f  = 0:1/t(end):.5/(t(2)-t(1));  %% positive frequency is f=0:1/T:0.5/(dt)
%%If = 2:41;


%If = 5:28;
If = 2:61;
If = 11:58;
If = 11:70;
%If = [11:34]+10;
%If = [11:34; %% 5-17Hz
% wavelet
model.f0 =10;
model.t0 =0;

% model
m  = 1e6./v(:).^2;
m0 = 1e6./v0(:).^2;


Q = speye(length(model.xsrc));
tic
% make data

R   = opRestriction(length(f),If);
%DFT = opKron(R*opDFTR(length(t)),opDirac(nd(2)),opDirac(nd(1)));
%DFT = opKron(R*opDFTR(length(t)),opDirac(nd(2)),opDirac(nd(1)));
%Dt =  DFT'*D;

%win = tukeywin(nd(3),0.4);
%D = vec(reshape(D,nd(1)*nd(2),nd(3))*diags(win));

%Dt = reshape(Dt,nd(1)*nd(2),nt);
disp('Preprocessing Filter is Completed......')


mask = zeros(n);
mask(2:end-1,5:end-6) = ones(n(1)-2,n(2)-10);
mask(151:end,:) = 0;
model.mask = mask(:);
model.mref = m;



%options.fid = fopen(['iter_awi_gaussian_lowhigh.log'],'w');

 options.itermax = 200;
% 
 options.tol = 1e-6;
 options.write = 1;
 options.maxit = 50;
 options.method ='lbfgs';

%% ouput the gather
disp('Generating Crosswell Data starts');
tic

model.freq = f(If);

D = F(m, Q, model);  D = gather(D);

toc
disp('Crosswell data is completed');

[od,dd,nd] = grid2odn(model.xrec,model.xsrc,model.freq);

R   = opRestriction(length(f),If);
DFT = opKron(R*opDFTR(length(t)),opDirac(nd(2)),opDirac(nd(1)));
%lambda = 1e12; alpha = 0.0;




ntrace = nd(1)*nd(2); nt = length(t); nfreq = length(model.freq);
Wd = tukeywin(nd(3),0.4); Wd = diags(Wd);
Wr = diags(tukeywin(nd(1),0.2));
Dt = DFT'*vec(Wr*reshape(reshape(D,ntrace,nd(3))*Wd,nd(1),nd(2)*nd(3))); 
Dt = reshape(Dt,prod(nd(1:2)),nt); figure,imagesc(1:ntrace,t,-Dt');


D0 = F(m0,Q,model); D0 = gather(D0);

Dp = D - D0;

Dpt = DFT'*vec(Wr*reshape(reshape(Dp,ntrace,nd(3))*Wd,nd(1),nd(2)*nd(3)));

Dpt = reshape(Dpt,prod(nd(1:2)),nt); figure,imagesc(1:ntrace,t,-Dpt');

idx = 2:(n(2)-1);
model.idx = idx; 


lambda = 1;  alpha = 0;
model.lambda = lambda; model.alpha = alpha;

%DisplayModel(model) 

model


opts.M = 5;
opts.maxit = 151;
opts.method = 'lbfgs';
opts.write = 1;
opts.vmin = 2.0; opts.vmax=5.5;
opts.tol = 1e-10;

    model.alpha = alpha;
    outinfo = DisplayModel(model);
    options.optTol = 1e-10; %% gradient;
    options.maxIter = 200;
    options.testOpt = 1;
    options.progTol = 1e-14;
    options.verbose  = 3;
    LB = 1/3.0^2*ones(size(m0)); UB=1/1.4^2*ones(size(m0));

%fprintf('Projected Gradient Method\n');
    fprintf('Projected LBFGS  Method\n');
   funProj = @(x)boundProject(x,LB,UB);
     tic
%curdir=['Gv5ReflectorsFWIfreq10to22Hz']
%mkdir(curdir);
%cd(curdir);
    %mf = minConf_SPG3(fh,x0,funProj,options);
    %mf = minConf_SPG2(fh,m0,funProj,options);
%   mf = minConf_PQN(fh,m0,funProj,options);

%mn = QGNewton(fh,m0,opts);
fh = @(x) MSWI_LineSrcfreqBatchTest(x, Dp,  model, Q, lambda, alpha); 
tic
[fval, F0, HH]=fh(m0);
F0 = gather(F0);
toc


R   = opRestriction(length(f),If);
nidx = length(idx); nt = length(t); 
nrec = nd(1); nsrc = nd(2); nfreq = nd(3); ntrace = nsrc*nrec;

DFT2 = opKron(R*opDFTR(length(t)),opDirac(nd(2)),opDirac(nidx));

F0 = F0*(Wd);
F0t = DFT2'*F0(:); F0t = reshape(F0t,nidx*nd(2),nt); F0t = fftshift(F0t,2);

figure, imagesc(1:ntrace,t-1,F0t');

% Ft0 =  datatransp(Ft,[nidx nd(2) nt],0);  Ft0 = reshape(Ft0,nidx,nd(2)*nt);
% figure, imagesc(Ft0);


Sft = DFT'*vec(HH.Sf*Wd); Sft = reshape(Sft,nd(1)*nd(2),length(t));
figure, imagesc(1:ntrace,t,-Sft');









