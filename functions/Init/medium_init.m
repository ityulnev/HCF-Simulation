%Create structure that contains all properties of the medium in which the
%pulse propagates, e.g. gas and fiber
classdef medium_init
    
   properties
   r_hcf,area_hcf,fullarea_hcf,Fluence,alpha,temperature,pressure,n2,n0,n,Iconst,gas,r_mode,k,k0,n0pressure,npressure,n_ext,npressure_ext,kph,kgr,kGV,kGVD,kTOD,k1_w0
   end
methods
    function s=medium_init(mesh,beam,gas)
    s.r_hcf=200e-6;%[m]  77.3345
    s.r_mode=0.65*s.r_hcf;%radius of zeroth mode *0.65
    s.area_hcf=pi*(s.r_mode)^2;%[m^2] effective mode area taken into account!  !!!as r=[0,R] take half area for calculations!!!
    s.fullarea_hcf=(s.r_hcf)^2*pi;
    % s.Fluence=beam.Q_In/s.area_hcf;
    s.alpha=log(beam.Q_In/beam.Q_Out);% 0.8
%     s.alpha=0;
    s.gas=gas;
    s.temperature=300;%[K]
    s.pressure=2;%[bar]

    switch gas
        case 'Neon'
        [s.n0pressure,s.n0]=calc_refrIndex(beam.wavelength,gas,s.pressure,s.temperature); %regular refractive index
        [s.npressure,s.n]=calc_refrIndex(mesh.wvl,gas,s.pressure,s.temperature); %regular refractive index
        [n2pressure,n2]=calc_refrIndex(beam.wavelength,'Neon_n2',s.pressure,s.temperature);
%         s.n2=2*7.5e-25;%@1barr pressure #### [m^2/W] Nonlinear refractive index
        s.n2=(n2pressure-s.n0pressure)/2.224e18;
        case 'Helium'
        s.n2=s.pressure*7.5e-25;%0.625e-24;%[m^2/W]          -> SPM  (work in progress value for neon)
        [s.n0pressure,s.n0]=calc_refrIndex(beam.wavelength,gas,s.pressure,s.temperature); %regular refractive index
        [s.npressure,s.n]=calc_refrIndex(mesh.wvl,gas,s.pressure,s.temperature); %regular refractive index
    end

    dfbound=find(mesh.f>beam.f0,1)-1-mesh.fbound;%shift f0 into the middle of array
    s.n_ext=[zeros(1,(mesh.fbound-1)-dfbound),s.n,zeros(1,dfbound)];
    s.npressure_ext=[zeros(1,(mesh.fbound-1)-dfbound),s.npressure,zeros(1,dfbound)];%[zeros(1,mesh.fbound-1),s.npressure];

    s.Iconst=0.5*s.n0*const.c*const.eps0;% Constant factor for calculating signal intensity I=Iconst*abs(E)^2
    % wave number k
    s.k0=s.n0pressure*beam.f0*2*pi/const.c;
    s.k=(s.npressure_ext).*(2.*pi.*mesh.f+beam.f0*2*pi)./const.c;
    % expanded k=k0+k'+k''+k'''  
    s.kGV=[0,diff(s.k)./(mesh.df*2*pi)]; 
    s.kGV((mesh.fbound)-dfbound)=0;
    s.kGV(mesh.flength-dfbound+1)=0;
    s.kGVD=[diff(s.kGV)./(mesh.df*2*pi),0];
    s.kGVD((mesh.fbound)-dfbound)=0;
    s.kGVD(mesh.flength-dfbound)=0;
    %         s.kTOD=[diff(s.kGVD)./(mesh.df*2*pi),0];
    s.k1_w0=s.kGVD(find(mesh.f>beam.f0,1)-1);
    
    end
end
end