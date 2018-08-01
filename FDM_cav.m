%% PhC L3 cavity simulation

RES=20; %resolution
dx=1/RES; dy=dx;
L=[20-dx/4 16*sqrt(3)/2-dy/4];% simulation Box size [Ly,Lx]
pml=[1 1 1 1]; % pml size on each boundary 
a=225; %Lattice parameter in nm
x2=-L(2):(dx/2):L(2);
y2=-L(1):(dx/2):L(1);
eps_h=1; % dielectric constant in the holes
lambda=1020/a; % wavelength of the simulation
eps_b=11.7;
r_a=61/a; % r/a radius of holes
f0=1/lambda; % central frequency for searching modes


type='triangular';
wgl=3;
h_r=zeros(1,wgl); % Array of holes sizes 
h_shx=zeros(1,wgl); % Array of holes x-shift along the array
h_shy=zeros(1,wgl); % Array of holes y-shift along the array
edge=[-1 -1 -1 -1];

averaging=11; % pixel averaging, should be an odd number
eps2=makephc(x2, y2, eps_h, eps_b, r_a, type, h_r, h_shx,h_shy,edge,averaging);%matrix of dielectric constant
mu2=ones(size(eps2));%matrix of permitivity
N_eig=10; % number of eigenvalues
figure
imagesc(eps2)
title('PhC cavity')

TE=true;
[WF f]=solve_fdm(dx,dy,f0,eps2,mu2,N_eig,TE);

%plot the Ey field of the M0 mode
figure
imagesc(abs(WF(:,:,2,6)))
title('Mode of an L3 cavity')


