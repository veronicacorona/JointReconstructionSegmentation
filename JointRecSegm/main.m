%% Joint reconstruction and segmentation from undersampled MRI data
% Veronica Corona, DAMTP, Cambridge
% vc324@cam.ac.uk
%% Load phantom image
load('brainphantom.mat')

%% Fourier operator
F = mridft(size(gt), [(size(gt,1)/2)+1 (size(gt,1)/2)+1]);

%% SPIRAL SAMPLING
load spiralsampling.mat %15%
S=speye(numel(samp));
S=S(samp(:),:);
nz = size(nonzeros(samp),1)/numel(samp); % undersampling rate
SF = S*F;
g = SF*gt(:);

% Gaussian noise
SIGMA = 0.005;
n = SIGMA*max(abs(g)).*randn(size(g));
f = g + n;

%% Gradient  operator
ex = ones(size(gt, 2), 1);
ey = ones(size(gt, 1), 1);
dx = spdiags([-ex ex], 0:1, size(gt, 2), size(gt, 2));
dx(size(gt, 2), :) = 0;
dy = spdiags([-ey ey], 0:1, size(gt, 1), size(gt, 1));
dy(size(gt, 1), :) = 0;

Grad = [kron(dx, speye(size(gt, 1))); kron(speye(size(gt, 2)), dy)];

%% Plot Zero-filling reconstructions
figure
%set(figure,'defaulttextinterpreter','latex');
subplot(1, 4, 1), imagesc(gt); axis off; title('Groundtruth')
subplot(1, 4, 2), imagesc(samp, [0 1]); axis off; title(['Undersampling matrix ', num2str(nz)])
subplot(1, 4, 3), imagesc(reshape(abs(SF'*g), size(gt)), [0 ...
    max(abs(gt(:)))]); axis off; title('Zero-filled reconstruction of noise-free data')
subplot(1, 4, 4), imagesc(reshape(abs(SF'*f), size(gt)), [0 ...
    max(abs(gt(:)))]); axis off ; title('Zero-filled reconstruction of noisy data')
colormap(gray)

%%  TV regularised Reconstruction
clear solver
alpha0=0.2;
K = alpha0*Grad;
sigma0=1/normest(K);
tau0=1/normest(K);

G = dataterm(S, F);
Fstar = projection(size(gt));
solver = pdhgm(K, Fstar, G);
G.setproxparam(tau0)
Fstar.setproxparam(sigma0)
solver.setmaxiter(1000)
solver.settolerance(5*10^-4)
solver.disableplot

G.setproxdata(f);
solver.solve;
u0=reshape(real(solver.getvariables.x), size(gt));
rel_tvrec = norm(gt - u0)/norm(gt);
figure;imagesc(u0);axis off;axis image;colormap(gray)
title(['TV Reconstruction, alpha=' num2str(alpha0) ', RRE =' num2str(rel_tvrec)])

%% Segmentation on TV regularised Reconstruction
beta0 = 0.001; %regulatisation parameter
c1 = 0.01; c2 = 0.3; c3 = 0.65; c4 = 0.8; % segmentation constants
classes = [c1 c2 c3 c4]; classes = classes';

vd1=convex_segmentation(u0, beta0, classes);

v02 = zeros(size(gt));
v02(vd1(:,:,1)==1)=1;
v02(vd1(:,:,2)==1)=2;
v02(vd1(:,:,3)==1)=3;
v02(vd1(:,:,4)==1)=4;

RSE = size(find(gt_seg~=v02),1)/numel(gt_seg);
figure; imagesc(v02); title(['Segmentation on TV reconstruction, \beta=' num2str(beta0) ', RSE =' num2str(RSE)])
%% Bregman TV Reconstruction
alpha01=1.1; %regularisation parameter
K01 = alpha01*Grad;
sigma0=1/normest(K01);tau0=1/normest(K01);
pk=zeros(size(gt));pk = pk(:);
figure
ulast=zeros(size(gt));
u01=ulast;
i=1;
while norm(SF*u01(:)-f(:)) > 0.005*max(abs(g)) * sqrt(numel(f))
    ulast=u01;
    
    G = dataterm_rec_bregman(S, F);
    Fstar = projection(size(gt));
    
    solver = pdhgm(K01, Fstar, G);
    G.setproxparam(tau0)
    Fstar.setproxparam(sigma0)
    solver.setmaxiter(1000)
    solver.settolerance(5*10^-4)
    solver.disableplot
    
    G.setproxdata(f);
    G.pk=pk;
    solver.solve;
    u01=reshape(real(solver.getvariables.x), size(gt));
    pklast=pk;
    pk=pk-(1/alpha01)*(real(F'*S'*(SF*u01(:)-f)));
    subplot(5,8,i); imagesc(u01);axis off;colormap gray
    i=i+1;
end
u0B=ulast;
RRE_breg=norm(gt - u0B)/norm(gt);
figure;imagesc(u0B);axis off;axis image;colormap(gray)
title(['Bregman Reconstruction, alpha=' num2str(alpha01) ', RRE =' num2str(RRE_breg)])
%% Segmentation on Bregman TV reconstruction
vB = convex_segmentation(u0B, beta0, classes);

v0B = zeros(size(gt));
v0B(vB(:,:,1)==1)=1;
v0B(vB(:,:,2)==1)=2;
v0B(vB(:,:,3)==1)=3;
v0B(vB(:,:,4)==1)=4;

RSE_breg = size(find(gt_seg~=v0B),1)/numel(gt_seg);
figure; imagesc(v0B); title(['Segmentation on Breg reconstruction, \beta=' num2str(beta0) ', RSE =' num2str(RSE_breg)])

%% Joint reconstruction and segmentation
%% Parameter initialization
alpha = 0.8; beta = 10e-7;
c1 = 0.01; c2 = 0.3; c3 = 0.7; c4 = 0.85;
classes = [c1 c2 c3 c4]; classes = classes';

delta = 0.001;

K1 = alpha*Grad;
K2 = beta*Grad;

sigma=1/normest(K1);tau=1/normest(K1);
sigma1=1/normest(K2);tau1=1/normest(K2);
y0 = 0.1;

v=v02(:);
V=0;
vd = vd1;
for i  = 1:length(classes)
    V=V+vd(:,:,i);
end

V1=0;
for i  = 1:length(classes)
    V1=V1+vd(:,:,i).* classes(i);
end
A =delta*V(:);
B = delta*V1(:);
rec(:,:,1)=u0; rec(:,:,2)=zeros(size(gt));
seg(:,:,1)=reshape(v,size(gt));
u = zeros(size(gt));
K1 = alpha*Grad;
sigma=1/normest(K1);
tau=1/normest(K1);
%% Joint reconstruction and segmentation
pk=zeros(size(gt)); pk = pk(:);
qk=zeros([numel(gt) length(classes)]);

figure
N=6;

for i=1:N
    clear solver1
    v=v(:);
    V=0;
    for m = 1:length(classes)
        V=V+vd(:,:,m);
    end
    
    V1=0;
    for m = 1:length(classes)
        V1=V1+vd(:,:,m).* classes(m);
    end
    A =delta*V(:);
    B = delta*V1(:);
    Gu = datatermBregman(S, F, A, B);
    Fu_star = projection(size(gt));
    Gu.setproxparam1(alpha);
    Gu.setproxdata(f);
    Gu.setproxparam(tau);
    Fu_star.setproxparam(sigma);
    Gu.pk=pk;
    clear solver1
    solver1 = pdhgm(K1, Fu_star, Gu);
    solver1.setmaxiter(1000);
    solver1.settolerance(5*10^-4);
    solver1.disableplot;
    solver1.sens=0.001;
    solver1.solve;
    
    u=reshape(real(solver1.getvariables.u),size(gt));
    pk=pk-(1/alpha)*(real(F'*S'*(SF*u(:)-f))-A.*u(:)-B);
    
    rec(:,:,i)=u;
    rel_rec_error(i) = norm(gt - u)/norm(gt);
    
    vlast=vd;
    delta1=delta;
    [~, ~, us, vd]=convex_segmentation_bregman(u, beta, classes, delta1, qk);
    v = zeros(size(gt));
    v(vd(:,:,1)==1)=1;v(vd(:,:,2)==1)=2;v(vd(:,:,3)==1)=3;v(vd(:,:,4)==1)=4;
    qk(:,1)=qk(:,1)-(1/beta)*delta1*((u(:)-c1)).^2;
    qk(:,2)=qk(:,2)-(1/beta)*delta1*((u(:)-c2)).^2;
    qk(:,3)=qk(:,3)-(1/beta)*delta1*((u(:)-c3)).^2;
    qk(:,4)=qk(:,4)-(1/beta)*delta1*((u(:)-c4)).^2;
    rel_seg_error(i) = size(find(gt_seg~=v),1)/numel(v);
    
    seg(:,:,i)=v;
    subplot(2,N,i); imagesc(u);title(['M= ' num2str(rel_rec_error(i))]); colormap gray
    subplot(2,N,i+N); imagesc(v); title(['M= ' num2str(rel_seg_error(i))]);
    
end
%% Plot results
figure
subplot(231)
imagesc(u0);axis off;axis image;colormap(gray)
title(['RRE=' num2str(rel_tvrec)])
subplot(232)
imagesc(u0B);axis off;axis image;colormap(gray)
title(['RRE=' num2str(RRE_breg)])
subplot(233)
imagesc(rec(:,:,N));title(['RRE= ' num2str(rel_rec_error(N))]); axis off;colormap gray
subplot(234)
imagesc(v02); title(['RSE=' num2str(RSE)]); axis off
subplot(235)
imagesc(v0B);axis off;axis image;colormap(gray)
title(['RSE=' num2str(RSE_breg)])
subplot(236)
imagesc(seg(:,:,N));title(['RSE= ' num2str(rel_seg_error(N))]);axis off;colormap gray