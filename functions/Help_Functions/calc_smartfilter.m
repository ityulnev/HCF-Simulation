%Function to build a smoothing filter on x which is:
%0:[1:x0-1]
%tanh(x):[x0:xmid-1]
%supergaussian(x):[xmid:end]
function [Sfilter]=calc_smartfilter(mesh,y1,pulse)
% xmid=find(max(abs(y1).^2)==(abs(y1).^2),1);
LRbounds=find_bounds(y1);
xright=LRbounds(1,2);
Gfilter=calc_supergaussian(mesh.f,10*mesh.df*(xright-pulse.pfmid),10,pulse.w0/(2*pi));
Tfilter=calc_tanhfilter((2*pi.*mesh.f-pulse.w0).*mesh.t,mesh.fmid);

Sfilter=[Tfilter(1:pulse.pfmid-1),Gfilter(pulse.pfmid:end)];
end