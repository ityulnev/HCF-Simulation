%Eliminate linear spectral phase components by shifting the max of electric
%field Et to t0=0 in time domain
%
function [Ef,t_delay]=cmpns_tshift(Ef,mesh)
Et=myifft(Ef,mesh);
[~,index_tmax]=max(abs(Et).^2); 
t_delay=(index_tmax-mesh.fmid).*mesh.dt;
Ef=Ef.*exp(1i*2*pi.*t_delay.*mesh.f);

end