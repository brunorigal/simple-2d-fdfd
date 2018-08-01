function [EM,f] = solve_fdm(dx,dy,f0,eps2,mu2,N_eig,TE,kinc)

% solve_fdm is an FDM EM solver on a Yee grid
% It returns WF, the Hz or Ez components of the wavefunctions and f their 
% frequency in inverse unit of length
% dx and dy are spatial resolution of the single grid (not the double grid)
% f0 is the initial guess for frequency
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

k0=2*pi*f0;
BC=[-2 -2]; % periodic boundary conditions
[Nx2, Ny2]=size(eps2);
[DEX,DEY,DHX,DHY] = diff_yee2(([Nx2 Ny2]/2),[dx dy],BC,kinc);

if TE % TE mode

epsxx=eps2(2:2:Nx2,1:2:Ny2);epsxx_inv=diag(sparse(1./epsxx(:)));
epsyy=eps2(1:2:Nx2,2:2:Ny2);epsyy_inv=diag(sparse(1./epsyy(:)));
muzz=mu2(2:2:Nx2,2:2:Ny2);muzz_inv=diag(sparse(1./muzz(:)));
H=-muzz_inv*(DEX*epsyy_inv*DHX+DEY*epsxx_inv*DHY);

else % TM mode
epszz=eps2(1:2:Nx2,1:2:Ny2);epszz_inv=diag(sparse(1./epszz(:)));
muxx=mu2(1:2:Nx2,2:2:Ny2);mu_x_inv=diag(sparse(1./muxx(:)));
muyy=mu2(2:2:Nx2,1:2:Ny2);mu_y_inv=diag(sparse(1./muyy(:)));
H=-epszz_inv*(DHX*mu_y_inv*DEX+DHY*mu_x_inv*DEY);
end
    
tic
[WFs,ks]=eigs(H,N_eig,k0^2); % solver
toc

[k, order] = sort(sqrt(diag(ks)));

eta_0=376.73031346177; % free space impedance
for i=1:length(ks)
    if TE
        EM(:,1,i)=WFs(:,order(i))*1i/eta_0; %Hz
        EM(:,2,i)=epsxx_inv*DHY*WFs(:,order(i))/k0; %Ex
        EM(:,3,i)=-epsyy_inv*DHX*WFs(:,order(i))/k0; %Ey
    else
        EM(:,1,i)=WFs(:,order(i)); %Ez
        EM(:,2,i)=muxx_inv*DHY*WFs(:,order(i))*1i/eta_0/k0; %Hx
        EM(:,3,i)=-muyy_inv*DHX*WFs(:,order(i))*1i/eta_0/k0; %Hy
    end
end
EM=reshape(EM,[Nx2/2, Ny2/2,3,N_eig]);
f=k/(2*pi);
end

