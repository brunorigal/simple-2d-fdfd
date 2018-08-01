function [Pix,Piy ] = Poynting_TE( EM,TE )
% returns the x and y components of the Poynting vector
% TE=true is the fdfd modes are TE

 if TE
    Pix(:,:)=real((EM(1:(end-1),:,3)+EM(2:end,:,3))/4.*conj(EM(1:(end-1),:,1)));
    Piy(:,:)=-real((EM(:,1:(end-1),2)+EM(:,2:end,2))/4.*conj(EM(:,1:(end-1),1)));
 else
    Pix(:,:)=-real(EM(:,:,3).*conj((EM(1:(end-1),:,1)+EM(2:end,:,1))/4));
    Piy(:,:)=real(EM(:,:,2).*conj((EM(:,1:(end-1),1)+EM(:,2:end,1))/4));
 end
end

