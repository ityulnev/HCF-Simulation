%Energy conserving fft  (t->f)
function [Ef]=myfft(Et,mesh)
if size(Et,1)>1
    Ef=fftshift(fft(ifftshift(Et),[],2))./(length(Et)*mesh.df);
    if strcmp(mesh.test,'yes')
        intEt=trapz(mesh.t,abs(Et(1,:)).^2);
        intEf=trapz(mesh.f,abs(Ef(1,:)).^2);
    end
elseif size(Et,1)==1
    Ef=fftshift(fft(ifftshift(Et)))./(length(Et)*mesh.df);
    if strcmp(mesh.test,'yes')
        intEt=trapz(mesh.t,abs(Et).^2);
        intEf=trapz(mesh.f,abs(Ef).^2);
    end
else
    error('myfft: size m of mxn matrix not real pos integer')
end
%Energy conservation check
if strcmp(mesh.test,'yes')
    tolerance=1e-4;
    test_errorMSG(abs(intEf-intEt)/intEf >tolerance,'myfft: Energy of Et and Ef not conserved!')
end  
end
