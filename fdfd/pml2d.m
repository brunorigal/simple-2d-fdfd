
function [sx,sy] = pml2d(NGRID,NPML,pml)

% pml2d computes the pml parameters sx and sy
% to absorb escaping waves on a 2D grid.
% NGRID= [ Nx Ny ] number of points in the grid
% NPML= [ Nx_low Nx_high Ny_low Ny_high ]size of the PML at the boundaries
% original parameters:
if nargin<3
    a_max=3;
    sig_max=1;
    p=3;
    type='norm';
elseif strcmp(pml.type,'PhC')
    type=pml.type;
    sig_max=pml.sig_max;
    a_max=0;
p=0;
else
    type='norm';
    a_max=pml.a_max;
    sig_max=pml.sig_max;
    p=pml.p;
end

% parameters optimised for PhC:
% a_max=10;
% sig_max=0.01; % maximum loss

sx=ones(NGRID);
sy=ones(NGRID);

sx(1:NPML(1),:)=kron(ones(1,NGRID(2)),(pml1d(NPML(1),a_max,sig_max,p,type).'));
sx((NGRID(1)-NPML(2)+1):NGRID(1),:)=flipud(kron(ones(1,NGRID(2)),pml1d(NPML(2),a_max,sig_max,p,type).'));
sy(:,1:NPML(3))=kron(ones(NGRID(1),1),(pml1d(NPML(3),a_max,sig_max,p,type)));
sy(:,(NGRID(2)-NPML(4)+1):NGRID(2))=kron(ones(NGRID(1),1),fliplr(pml1d(NPML(4),a_max,sig_max,p,type)));
end


function [s] = pml1d(NPML,a_max,sig_max,p,type)
eta_0=376.73031346177; % free space impedance
if strcmp(type,'PhC');
    s=1+1i*eta_0*sig_max*exp(1-1./linspace(1,1/NPML,NPML));
else 
    s=(1+a_max*linspace(1,1/NPML,NPML).^p).*(1+1i*eta_0*sig_max*sin(pi/2*linspace(1,1/NPML,NPML)).^2);
end
end
