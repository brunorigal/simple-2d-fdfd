function [eps] =makephc(x, y, eps_h, eps_b, r_a, type, h_r, h_shx,h_shy,edge,averaging)
% the output eps is a matrix with the refractive index at each pixel
% x and y are the x and y coordinates of the simulation.
% eps_h is the refractive index of holes.
% eps_b is the refractive index of the bulk material.
% r_a is r/a the radius of not changed holes in units of a. 
% h_r, h_sh1, h_sh2 are matrixes containing the holes radia and shifts in
% units of a.
% edge is the distance to the edge of the simulation without holes.

dx=x(2)-x(1);
dy=y(2)-y(1);
    if nargin>10
        av=averaging; % number of pixels of the averaging subgrid per pixel. It should be an odd number.
        lav=floor(av/2);
        xf=linspace((min(x)-0.5*dx),(max(x)+0.5*dx),(length(x)+1)*av); 
        yf=linspace((min(y)-0.5*dy),(max(y)+0.5*dy),(length(y)+1)*av); 
    end
if strcmp(type, 'triangular')
    a1 = [0 1];
    a2 = 0.5*[sqrt(3) 1]; 
elseif strcmp(type, 'square')
    a1 = [0 1];
    a2 = [1 0]; 
end

[Y0,X0] = meshgrid(y,x);

    eps = eps_b .* ones(size(X0)); 
    [SX,SY]  = size(eps);

[Smx,Smy]=size(h_r);

for j=-80:80
    for i=-80:80
        posj=(j+floor(Smx/2)+1);
        posi=(i+floor(Smy/2)+1);
         if (posj<=Smx)&&(posj>0)&&(posi<=Smy)&&(posi>0)
                XA = a1(1)*i + a2(1)*j+h_shx(posj,posi);                
                YB = a1(2)*i + a2(2)*j-h_shy(posj,posi);
                ra=h_r(posj,posi)*r_a;
        else
                XA = a1(1)*i + a2(1)*j;
                YB = a1(2)*i + a2(2)*j;
                ra=r_a;
        end
        if ((XA>min(x+edge(1)))&(XA<max(x-edge(2)))&(YB>min(y+edge(3)))&(YB<max(y-edge(4))))

            R = sqrt((X0-XA).^2+(Y0-YB).^2);
            if nargin>10
                % for scalar eps averaging
                [xh,yh]=find(R==min(min(R)),1); 
                hs=floor(ra/dx)+2; 
                ys=max((-hs+yh),1):min((hs+yh),length(y));
                xs=max((-hs+xh),1):min((hs+xh),length(x)); 
                [Y0s,X0s] = meshgrid( ys,xs);

                xsf=(min(xs)*av-lav):(max(xs)*av+lav);
                ysf=(min(ys)*av-lav):(max(ys)*av+lav);
                [Y0f,X0f] = meshgrid(yf(ysf),xf(xsf));
                Rf= sqrt((X0f-XA).^2+(Y0f-YB).^2);
                mark= ones(length(xsf),length(ysf));
                mark(find(Rf<ra))=0;
                
                mark_av=squeeze(sum(sum(reshape(mark,av,size(mark,1)/av,av,size(mark,2)/av),1),3)/av^2);
                eps(xs,ys)=eps_h+(eps_b-eps_h)*mark_av;                
                                
            else
                eps(find(R<ra)) = eps_h;
            end  
     end
        end
    end
end