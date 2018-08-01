%% PhC waveguide simulation
close all; clear all;
res=20; % resolution
dx=1/res; dy=dx;
L=[0.5-dx/4 16*sqrt(3)/2+dy/4 ];% simulation Box size [Ly,Lx]
pml=[1 1 1 1]; %pml size on each boundary 
a=225; %Lattice parameter in nm
x2=-L(2):(dx/2):L(2);
y2=-L(1):(dx/2):L(1);

lambda=1000/a;
eps_h=1; % dielectric constant in the holes
eps_b=11; % dielectric constant of the material: GaAs
r_a=61/a; % r/a radius of holes
f0=1/lambda; % central frequency for searching modes


type='triangular';
wgl=7;
h_r=zeros(1,wgl); % Array of holes sizes 
h_shx=zeros(1,wgl); % Array of holes x-shift along the array
h_shy=zeros(1,wgl); % Array of holes y-shift along the array

edge=[-1 -1 -1 -1];
averaging=11;

eps2=makephc(x2, y2, eps_h, eps_b, r_a, type, h_r, h_shx,h_shy,edge);
figure
imagesc(eps2)
title('PhC waveguide')

mu2=ones(size(eps2));
N_eig=20;
kvec=2*pi*(0:0.05:0.5);
for j=1:length(kvec)
    kinc=[kvec(j) 0];
    TE=true;
    [WF f(:,j)]=solve_fdm(dx,dy,f0,eps2,mu2,N_eig,TE,kinc);
end

figure
plot(kvec/(2*pi),real(f),'ob')
title('Band structure')

figure
imagesc(abs(WF(:,:,2,19)))
daspect([1,1,1])
title('Ey component of the waveguide mode')


