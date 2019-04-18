%Energy conserving ifft (f->t)
function [Et]=myifft(Ef,mesh)
if size(Ef,1)>1
    Et=fftshift(ifft(ifftshift(Ef),[],2)).*(length(Ef)*mesh.df);
elseif size(Ef,1)==1
    Et=fftshift(ifft(ifftshift(Ef))).*(length(Ef)*mesh.df);  
else
    error('myifft: Dimension m of mxn matrix not real pos integer')
end


If=abs(Ef).^2;
It=abs(Et).^2;
intEf=trapz(mesh.f,If);
intEt=trapz(mesh.t,It);
        %Energy conservation check
        tolerance=1e-4;
        if abs(intEf-intEt)/intEf <tolerance
            %do nothing  
        else
            warning('pulse_init: Energy of Et and Ef not conserved!')    
        end
          
end
