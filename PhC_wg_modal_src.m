%% Modal or dipole source in a PhC waveguide

% Make the geometry 
% Warning!!! the vertical direction is y and the horizontal direction is x!
res=20; % simulation resolution
dx=1/res; dy=dx;
L=[20 16*sqrt(3)/2]; % simulation Box size [Ly,Lx]
L=2*floor(L/(dx/2)/2)*dx/2+dx/4;
pml=[5 5 1 1 ]; % pml length [l_y_left,l_y_right,l_x_down,l_x_up] on each boundary 
BC=[-2 -2]*0; % 0->Dirichlet Boundary conditions
TE=true; % Transverse electric simulation

x=-L(2):dx:L(2);
y=-L(1):dy:L(1);
x2=-L(2):(dx/2):L(2);
y2=-L(1):(dy/2):L(1);
eps_h=1; % dielectric in the air
a=225; % lattice constant in nm.
r_a=61/a; % hole radius
lam_nref=1060/a;% wavelength used for the effective index
eps_b=11.6; % dielectric in the PhC material

type='triangular';
wgl=50; 
h_r=zeros(1,wgl); % Array of holes sizes 
h_shx=zeros(1,wgl); % Array of holes x-shift along the array
h_shy=zeros(1,wgl); % Array of holes y-shift along the array
edge=[1 1 1 1]*-0.1;  % edge with no holes
averaging=11;
eps2=makephc(x2, y2, eps_h, eps_b, r_a, type, h_r, h_shx,h_shy,edge,averaging);
figure
imagesc(eps2)
title('PhC waveguide')
mu2=ones(size(eps2));

% Modal source
source.type='modal_source_x';
f0=1/lam_nref; %search modes around this wavelength
kinc=2*pi*[0.4 0]; % wavevector of the input mode
posx=-8; % position of the input mode
pos(1)=findc(posx,x);
pos(2)=pos(1)+res-1;
N_eig=5; % number of modes to be searched
[WFmodal fmodal]=solve_fdm(dx,dy,f0,eps2(:,(2*pos(1)-1):(2*pos(2))),mu2(:,(2*pos(1)-1):(2*pos(2))),N_eig,TE,kinc);
nmode=find(fmodal>0.21,1); % select only the waveguide mode
source.Hz=zeros(size(eps2)/2); 
source.Hz(:,pos(1):pos(2))=WFmodal(:,:,1,nmode);
source.pos=pos; % position of the modal source
k0=2*pi*fmodal(nmode); % wavevector of the simulation
%
% % Hz dipole source
% source.type='point_source';
% source.pos=[0 0];
% source.radius=0.2; % radius of the masking matrix
% k0=2*pi/lam_nref;

print=1; % 1-> print only the resulting fields
[EM tim]=solve_fdfd(x,y,k0,eps2,mu2,pml,BC,source ,TE,print);

[Pix Piy]=Poynting_TE(EM,TE); % output the Poynting vectors
Poynt_y=sum(Piy,1)*dx;  
figure; plot(y(2:end),Poynt_y);
R=abs(Poynt_y(findc(-5,y))/Poynt_y(findc(5,y))); % reflection ratio
