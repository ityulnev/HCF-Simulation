%Energy conserving ifft (t->f)
function [Ef]=myifft(Et,dt)
if size(Et,1)>1
    Ef=fftshift(ifft(ifftshift(Et),[],2)).*(length(Et)*dt);
elseif size(Et,1)==1
    Ef=fftshift(ifft(ifftshift(Et))).*(length(Et)*dt);  
else
    error('myifft: Dimension m of mxn matrix not real pos integer')
end

end
