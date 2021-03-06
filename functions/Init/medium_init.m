%Create structure that contains all properties of the medium in which the
%puse propagates, e.g. gas and fiber
classdef medium_init
    
   properties
   r_hcf,area_hcf,fullarea_hcf,Fluence,alpha,temperature,pressure,
   n2,Beta2,n0,n,Iconst,gas,r_mode,k,k0,n0pressure,npressure,n_ext,npressure_ext,kph,kgr,kGV,kGVD,kTOD
   end
methods
    function s=medium_init(mesh,beam,gas)
    s.r_hcf=200e-6;%[m]
    s.r_mode=0.65*s.r_hcf;%radius of zeroth mode
    s.area_hcf=pi*(s.r_mode)^2;%[m^2] effective mode area taken into account!
    s.fullarea_hcf=(s.r_hcf)^2*pi;
    % s.Fluence=beam.Q_In/s.area_hcf;
    s.alpha=log(beam.Q_In/beam.Q_Out);% 0.8
%     s.alpha=0;
    s.gas=gas;
    s.temperature=300;%[K]
    s.pressure=2;%[bar]

    switch gas
        case 'Neon'
        s.n2=s.pressure*0.85e-24;%0.625e-24;%[m^2/W]                   -> SPM
        s.Beta2=-s.pressure*3.7*(1e-30);%[s^2/m]             -> GVD
        [s.n0pressure,s.n0]=calc_refrIndex(beam.wavelength,gas,s.pressure,s.temperature); %regular refractive index
        [s.npressure,s.n]=calc_refrIndex(mesh.wvl,gas,s.pressure,s.temperature); %regular refractive index
    end

    s.n_ext=[zeros(1,mesh.fbound-1),s.n];
    s.npressure_ext=[zeros(1,mesh.fbound-1),s.npressure];

    s.Iconst=0.5*s.n0*const.c*const.eps0;%Constant factor for calculating signal intensity I=Iconst*abs(E)^2
    s.k=s.npressure_ext.*(2.*pi.*mesh.f)./const.c;

    s.k0=s.n0pressure*beam.f0*2*pi/const.c;
    %% k 
%     dndw=[diff(s.n)./(mesh.df*2*pi),0];
%     dndw_ext=[zeros(1,mesh.fbound-1),dndw];
%     s.kph=s.n_ext./const.c;
%     s.kgr=s.kph+dndw_ext.*((2*pi.*mesh.f)./const.c);
    s.kGV=[diff(s.k)./(mesh.df*2*pi),0];
    s.kGVD=[diff(s.kGV)./(mesh.df*2*pi),0];
    s.kTOD=[diff(s.kGVD)./(mesh.df*2*pi),0];
    
    end
end
end