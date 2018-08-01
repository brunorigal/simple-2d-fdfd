function [WF,n0,kap] = solve_fdm_wg(dx,dy,nguess,lambda,eps2,mu2,N_eig,kinc)

% solve_fdm_wg is an FDM EM solver on a Yee grid for a 3D waveguide, with
% waves propagating along z.
% It returns WF, the Hz or Ez components of the wavefunctions and f their 
% frequency in inverse unit of length
% dx and dy are spatial resolution of the single grid (not the double grid)
% nguess is the initial guess for the modes effective index
% lambda is the wavelength
% eps2 is the dielectic on a double grid, the real simulated grid has 
% two times less resolution
% mu2 is the permeability on a double grid, the real simulated grid has 
% two times less resolution
% N_eig is the number of eigenvalues searched
% TE=true if the mode is TE and false if the mode is TM
% kinc is the wavevector associated with Bloch mode boundary conditions 
% (Periodic boundary conditions with a phase term)

if nargin<8
    kinc=[0 0];
end

k0=2*pi/lambda;
BC=[-2 -2 -2 -2]*0; % 0->Dirichlet Boundary conditions
[Nx2, Ny2]=size(eps2);

%Build the differential operators
[DEX,DEY,DHX,DHY] = diff_yee2(([Nx2 Ny2]/2),k0*[dx dy],BC,kinc);
epsxx=eps2(2:2:Nx2,1:2:Ny2);epsxx=diag(sparse(epsxx(:)));
epsyy=eps2(1:2:Nx2,2:2:Ny2);epsyy=diag(sparse(epsyy(:)));
epszz=eps2(1:2:Nx2,1:2:Ny2);epszz_inv=diag(sparse(1./epszz(:)));
muxx=mu2(1:2:Nx2,2:2:Ny2);muxx=diag(sparse(muxx(:)));
muyy=mu2(2:2:Nx2,1:2:Ny2);muyy=diag(sparse(muyy(:)));
muzz=mu2(2:2:Nx2,2:2:Ny2);muzz_inv=diag(sparse(1./muzz(:)));

P=[DEX*epszz_inv*DHY -(DEX*epszz_inv*DHX+muyy);
    DEY*epszz_inv*DHY+muxx -DEY*epszz_inv*DHX];
Q=[DHX*muzz_inv*DEY -(DHX*muzz_inv*DEX+epsyy);
    DHY*muzz_inv*DEY+epsxx -DHY*muzz_inv*DEX];
H=P*Q;

    
tic
[WFs,ks]=eigs(H,N_eig,-nguess^2); % solver
toc

[k, order] = sort((diag(ks)));
WFs=reshape(WFs,Nx2/2, Ny2/2,2,N_eig);

for i=1:length(ks)
    WF(:,:,:,i)=WFs(:,:,:,order(i));
end

n0=imag(sqrt((k)));
kap=real(sqrt((k)));

end

