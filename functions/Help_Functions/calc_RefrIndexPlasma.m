%Calculates the part of refractive index corresponding to plasma
%(valid for higher ionization % of the gas)
%Lorentz-Lorentz Equation for mixture of gas+plasma
function [dn_p,nL]=calc_RefrIndexPlasma(mesh,beam,medium,n_e,n_gas)

n0=medium.n0pressure;
wp2=const.e^2.*n_e./(const.eps0*const.m_e);
np2=(1-wp2./(beam.w0.^2));
DensityFrac=n_e./n_gas;

C=DensityFrac.*(np2-1)./(np2+2)+(1-DensityFrac).*(n0^2-1)./(n0^2+2);

nL2=(1+2.*C)./(1-C);
nL=sqrt(nL2);
dn_p=-nL+n0;

end