%Energy conserving fft  (t->f)
function [Ef]=myfft(Et,mesh)
if size(Et,1)>1
    Ef=fftshift(fft(ifftshift(Et),[],2))./(length(Et)*mesh.df);
elseif size(Et,1)==1
    Ef=fftshift(fft(ifftshift(Et)))./(length(Et)*mesh.df);
else
    error('myfft: size m of mxn matrix not real pos integer')
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
