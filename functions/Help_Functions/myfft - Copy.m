%Energy conserving fft  (f->t)
function [Et]=myfft(Ef,dt)
if size(Ef,1)>1
    Et=fftshift(fft(ifftshift(Ef),[],2))./(length(Ef)*dt);
elseif size(Ef,1)==1
    Et=fftshift(fft(ifftshift(Ef)))./(length(Ef)*dt);
else
    error('myfft: size m of mxn matrix not real pos integer')
end

end
