function [EM,t] = solve_fdfd(x,y,k0,eps2,mu2,pml,BC,source,TE ,print,kinc)

% solve_fdfd is a FDFD electromagnetic solver on a Yee grid
% it is based on: Rumpf, Raymond C. ?Simple Implementation of Arbitrarily 
% Shaped Total-Field/Scattered-Field Regions in Finite-Difference 
% Frequency-Domain.? Progress In Electromagnetics Research B 36 (2012): 221?248.

% 1) It returns EM the EM field components results of the simulation
% for a TE simulation, EM(:,:,1) is Hz, EM(:,:,2) is Ex, EM(:,:,3) is Ey 
% for a TM simulation, EM(:,:,1) is Ez, EM(:,:,2) is Hx, EM(:,:,3) is Hy
% 2) dx and dy are spatial resolution of the single grid (not the double grid)
% 3) f0 is the source frequency in inverse unit of length
% 4) eps2 is the dielectic on a double grid, the real simulated grid has 
% two times less resolution
% 5) mu2 is the permeability on a double grid, the real simulated grid has 
% two times less resolution
% 6) pml=[pmlx_low pmlx_high pmly_low pmly_high] contains the pml size at each
% boundary
% 7) BC is the Boundary conditions 0: Dirichlet Boundary conditions, -2:
% Periodic Boundary conditions
% 8) Source contains the info about the source, it is a struct and can be:
% -> Source.type='point_source', in which case it also needs the fields:
% pos: the position of the source
% radius: the radius of the masking matrix
% -> Source.type='plane_wave', in which case it also needs the fields:
% 
% 8) TE=true if the simulation is TE and false if the simulation is TM
% 9) print=1 prints the sim result, print=2 prints also intermediate figures 
% shape of the source, dielectric grid...)
% 10) kinc is the wavevector associated with Bloch mode boundary conditions 
% (Periodic boundary conditions with a phase term)

if nargin<11
    kinc=[0 k0];
end
dx=x(2)-x(1);dy=y(2)-y(1);
[Nx2, Ny2]=size(eps2); N=[Nx2 Ny2]/2;

%make the PML
if isstruct(pml)
     NPML=floor([(pml.s(3:4))/dy,(pml.s(1:2))/dx]);NPML2=2*NPML;
    [sx,sy] = pml2d([Nx2 Ny2],NPML2,pml);
else
    NPML=floor([(pml(3:4))/dy,(pml(1:2))/dx]);NPML2=2*NPML;
    [sx,sy] = pml2d([Nx2 Ny2],NPML2);
end


if print==2; figure; imagesc(y,x,sqrt(abs(sx).^2+abs(sy).^2));daspect([1 1 1]);title pml; end;

% Make the differential operator
[DEX,DEY,DHX,DHY] = diff_yee2(N,k0*[dx dy],BC,kinc/k0);
if TE
eps2xx=sy.*eps2./sx; epsxx=eps2xx(2:2:Nx2,1:2:Ny2);epsxx_inv=diag(sparse(1./epsxx(:)));
eps2yy=sx.*eps2./sy; epsyy=eps2yy(1:2:Nx2,2:2:Ny2);epsyy_inv=diag(sparse(1./epsyy(:)));
mu2zz=sx.*sy.*mu2;muzz=mu2zz(2:2:Nx2,2:2:Ny2);muzz=diag(sparse(muzz(:)));
A=DEX*epsyy_inv*DHX+DEY*epsxx_inv*DHY+muzz; % TE mode
else
eps2zz=sx.*eps2.*sy; epszz=eps2zz(1:2:Nx2,1:2:Ny2);epszz=diag(sparse(epszz(:)));
mu2xx=mu2.*sy./sx;muxx=mu2xx(1:2:Nx2,2:2:Ny2);muxx_inv=diag(sparse(1./muxx(:)));
mu2yy=mu2.*sx./sy;muyy=mu2yy(2:2:Nx2,1:2:Ny2);muyy_inv=diag(sparse(1./muyy(:)));
A=DHX*muyy_inv*DEX+DHY*muxx_inv*DEY+epszz; % TM mode
end

% Make the source
if strcmp(source.type, 'plane_wave_y')
    kinc=source.kinc;
    [Y,X] = meshgrid(y,x);
    fsrc = exp(1i*(kinc(2)*X + kinc(1)*Y));
    if print==2
        figure; imagesc(y,x,real(fsrc));daspect([1 1 1]);title source
        figure; imagesc(y,x,real(eps2(1:2:end,1:2:end)));daspect([1 1 1]); title dielectric
    end
    fsrc=fsrc(:);
    % CONSTRUCT Q the masking matrix
    lim=source.lim;
    Q = Y<lim; 

    if print==2; figure; imagesc(y,x,Q);daspect([1 1 1]);  title('mask'); end;

    Q = diag(sparse(Q(:)));
    b=(Q*A-A*Q)*fsrc;

elseif strcmp(source.type, 'point_source')
    [Y,X] = meshgrid(y,x);
    pos=source.pos;
    radius=source.radius;
    R = sqrt((X-pos(2)).^2 + (Y-pos(1)).^2); 
    [nx, ny]=find(min(R(:))==R);
    n_s=sqrt(eps2(2*nx,2*ny)); 
    fsrc = exp(1i*k0*n_s*R)./sqrt(R); 
    if print==2
        figure; imagesc(y,x,real(fsrc));daspect([1 1 1]); title('source');
        figure; imagesc(y,x,real(eps2(1:2:end,1:2:end)));daspect([1 1 1]); title('dielectric');
    end
    fsrc=fsrc(:);
    % CONSTRUCT Q the masking matrix
    Q = R <radius; %Q(:,[1:(NPML(3)+2) (end-(NPML(4)+2)):end])=1; Q([1:(NPML(1)+2) (end-(NPML(2)+2)):end],:)=1;
        if print==2 figure; imagesc(y,x,Q);daspect([1 1 1]); title mask; end
    Q = diag(sparse(Q(:)));
    b=(Q*A-A*Q)*fsrc; 
elseif strcmp(source.type, 'point_source_Ey')
    [Y,X] = meshgrid(y,x);
    pos=source.pos;
    radius=source.radius;
    R = sqrt((X-pos(2)).^2 + (Y-pos(1)).^2); 
    [nx, ny]=find(min(R(:))==R);
    n_s=sqrt(eps2(2*nx,2*ny)); 
    fsrc = exp(1i*k0*n_s*R).*(Y-pos(1))./(R); 
    if print==2
        figure; imagesc(y,x,real(fsrc));daspect([1 1 1]); title('source');
        figure; imagesc(y,x,real(eps2(1:2:end,1:2:end)));daspect([1 1 1]); title('dielectric');
    end
    fsrc=fsrc(:);
    % CONSTRUCT Q the masking matrix
    Q = R <radius; 
        if print==2 figure; imagesc(y,x,Q);daspect([1 1 1]); title mask; end
    Q = diag(sparse(Q(:)));
    b=(Q*A-A*Q)*fsrc; 
    
elseif strcmp(source.type, 'modal_source_x')
    [Y,X] = meshgrid(y,x);
    pos=source.pos;
    fsrc =source.Hz;
    if print==2
        figure; imagesc(y,x,real(fsrc));daspect([1 1 1]); title('source');
        figure; imagesc(y,x,real(eps2(1:2:end,1:2:end)));daspect([1 1 1]); title('dielectric');
    end
    fsrc=fsrc(:);
    % CONSTRUCT Q the masking matrix
    Q = Y <y(floor((pos(1)+pos(2))/2)); %Q(:,[1:(NPML(3)+2) (end-(NPML(4)+2)):end])=1; Q([1:(NPML(1)+2) (end-(NPML(2)+2)):end],:)=1;
    if print==2;  figure; imagesc(y,x,Q);daspect([1 1 1]); title mask; end
    Q = diag(sparse(Q(:)));
    b=(Q*A-A*Q)*fsrc; 
end
  
% solving FDFD
tic;
f=A\b;
t=toc

eta_0=376.73031346177; % free space impedance

% shape the output E field
if TE
    EM(:,1)=f*1i/eta_0; %Hz
    EM(:,2)=epsxx_inv*DHY*f/k0; %Ex
    EM(:,3)=-epsyy_inv*DHX*f/k0; %Ey
else
    EM(:,1)=f; %Ez
    EM(:,2)=muxx_inv*DHY*f*1i/eta_0/k0; %Hx
    EM(:,3)=-muyy_inv*DHX*f*1i/eta_0/k0; %Hy
end
EM=full(reshape(EM,[N,3]));

%print the result
if print==2||print==1
if TE; tt={'H_z','E_x','E_y'};
else  tt={'E_z','H_x','H_y'}; end;
h=figure;
subplot(1,3,1);
imagesc(y,x,log(abs(EM(:,:,1)))); daspect([1 1 1]); title(tt{1});
subplot(1,3,2);
imagesc(y,x,real(EM(:,:,2))); daspect([1 1 1]); title(tt{2});
subplot(1,3,3);
imagesc(y,x,real(EM(:,:,3))); daspect([1 1 1]); title(tt{3});
set(gcf,'Position',[100 100 700 300])
end


