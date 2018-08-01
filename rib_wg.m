% Rib waveguide simulation
close all; clear all
RES=30; %resolution
dx=1/RES; dy=dx;
L=[3-dx/4 3-dy/4];% simulation Box size [Ly,Lx] in um

x2=-L(2):(dx/2):L(2);
y2=-L(1):(dx/2):L(1);

ncore=1.9;
ndown=1.52;
nup=1;
nguess=1.9; % index at which the modes are searched. Should be slightly above ncore

eps2=ones(length(x2),length(y2))*nup;
eps2(findc(0,x2):findc(0.25,x2),:)=ncore^2;
eps2(findc(0,x2):findc(0.85,x2),findc(-1,x2):findc(1,x2))=ncore^2;
eps2(1:findc(0,x2),:)=ndown^2;
mu2=ones(size(eps2));
N_eig=10; % number of eigenvalues
figure; imagesc(x2,y2,eps2);
title('rib waveguide')
lambda=0.9; % simulated wavelength in um
TE=true;

[WF,n0,kap]=solve_fdm_wg(dx,dy,nguess,lambda,eps2,mu2,N_eig,TE);

i=1
figure; subplot(1,2,1);imagesc(abs(WF(:,:,2,i)));daspect([1 1 1]);
title('mode of a rib waveguide')

