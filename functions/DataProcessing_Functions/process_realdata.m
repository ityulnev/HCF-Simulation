%Process data
%Brings data into same format(unit,even-spacing,length) as in simulation
function [varargout]=process_realdata(varargin)
f=varargin{end};
for m=1:(length(varargin)-1)/2
%Convert from wavelength to energy
[Ix,Iy]=convert_data(varargin{m+(m-1)},varargin{m+m},'forward');
%As there are no negative Intensities!
Iy(Iy<0)=0;
%Interpolate for evenly spaced Signal E,I(E)
[Ix2,Iy2]=even_data(Ix,Iy);
%Unit coversion of x-axis  [J->1/s]
Ix2e=Ix2./const.h;
%Extending the freq. domain so that Simulation doesn't hit boundaries
[Ix3,Iy3]=extendToZero(Ix2e,Iy2);
%match my meshed f domain
Iy4 = interp1(Ix3,Iy3,f);
Iy4(isnan(Iy4))=0;
%nomalize
norm_factor=1/abs(trapz(f,Iy4));
varargout{m+(m-1)}=f;
varargout{m+m}=Iy4.*norm_factor;
end

end