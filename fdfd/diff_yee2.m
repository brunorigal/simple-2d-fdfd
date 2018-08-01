function [DEX,DEY,DHX,DHY] = diff_yee2(NGRID,RES,BC,kinc)

% Input Arguments
% =================
% NGRID [Nx Ny] grid size
% RES [dx dy] grid resolution of the 1X grid
% BC [xbc ybc] boundary conditions
% -2: periodic (requires kinc)
%  0: Dirichlet
% kinc % [kx ky] incident wave vector
% This argument is only needed for periodic boundaries.
% Position on the grid: m = (ny - 1)*Nx + nx

if nargin==4
    kinc0=kinc;
end
dx=RES(1); dy=RES(2);
NX=NGRID(1); NY=NGRID(2); N=NX*NY;
Lamx=NX*dx;
Lamy=NY*dy;
if NX==1
    if BC(1)==0;
        DEX=sparse(N,N);
    else
    DEX=spdiags(1i*kinc0(2)*ones(N,1),0,N,N)
    end
else
    DEX=spdiags([-1./dx.*ones(N,1),1./dx.*ones(N,1)],[0 1],N,N);
    for ny=1:(NY-1)
        DEX(p(NX,ny,NX),p(NX+1,ny,NX))=0;  % Dirichlet BC 
    end
    if BC(1)==-2
        for ny=1:NY
            DEX(p(NX,ny,NX),p(1,ny,NX))=exp(1i*kinc0(2)*Lamx)/dx; % periodic BC
        end

    end
end

if NY==1
    if BC(2)==0;
        DEY=sparse(N,N);
    else
        DEY=spdiags(1i*kinc0(1)*ones(N,1),0,N,N)
    end
else
    DEY=spdiags([-1./dy.*ones(N,1),1./dy.*ones(N,1)],[0 NX],N,N);
    if BC(2)==-2
        for nx=1:NX
            DEY(p(nx,NY,NX),p(nx,1,NX))=exp(-1i*kinc0(1)*Lamy)/dy; % periodic BC
        end
    end
end

  DHX=-DEX';
  DHY=-DEY';
end
function [n] = p(nx,ny,NX)
     n=nx+NX*(ny-1);
end