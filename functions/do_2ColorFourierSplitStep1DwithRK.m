%Calculate electric field propagation with steps in 0:dz:Lz
%d/dz(E) from SVEA in 1D
%Implemented Effects: GVD, SPM+SST, XPM, Attenuation, Ionization
function [EtF,EtS]=do_2ColorFourierSplitStep1DwithRK(mesh,beam,medium,pulse1,pulse2)
h=mesh.dz;

index2f0=find(mesh.f>=(pulse2.w0/(2*pi)),1)-1;
% Lwalkoff=pulse2.t_fwhm/abs(medium.kGV(indexf0)-medium.kGV(index2f0)); % Walkoff length for 2 pulses copropagating
Twalkoff=abs(medium.kGV(pulse1.pfmid)-medium.kGV(index2f0)).*mesh.Lz;
EfF=abs(pulse1.Erf);% Fundamental pulse
EfS=pulse2.Erf;% Second Harmonic pulse

[n_e,Eg,n_gas]=calc_eDensityADK((pulse1.Ert),mesh,medium,pulse1,beam);
IonizLvl=max(n_e)/n_gas;

    for m=1:(mesh.zlength)
        %% Group Velocity difference        
        dGV=-1i.*(2*pi*mesh.f).*(medium.kGV(index2f0)-medium.kGV(pulse1.pfmid));
        dGV2=-1i.*(2*pi*mesh.f).*(1/const.c-medium.kGV(pulse1.pfmid));
        %% Group velocity dispersion GVD of Fundamental
        %Fundamental
        GVD=-1i.*(2*pi.*mesh.f-pulse1.w0).^2.*medium.kGVD(pulse1.pfmid)./2;
        GVD(1:mesh.fmid)=0;
        EfF=EfF.*exp((GVD).*h).*exp(-medium.alpha/2*h);%alpha -> Energy loss      
        EfF(isinf(EfF))=0;
    %             EfF(1:mesh.fmid)=0;
        %% Group velocity dispersion GVD of Second Harmonic
        GVD2=-1i.*(2*pi.*mesh.f-pulse2.w0).^2.*medium.kGVD(index2f0)./2;
        EfS=EfS.*exp((GVD2+dGV).*h).*exp(-medium.alpha/2*h);%alpha -> Energy loss                     
        EfS(isinf(EfS))=0;
        %% Runge Kutta with Self phase modulation SPM and XPM ...
        newEfF=do_oneRKstep(mesh,pulse1,pulse2,beam,medium,EfF,EfS,m,h);
%         newEfS=do_oneRKstep(mesh,pulse2,pulse1,beam,medium,EfS,EfF,m,h);
        EfF=newEfF;
%         EfS=newEfS;        
        EtF=myifft(EfF,mesh); 
        EtS=myifft(EfS,mesh);   
            %% Plot
            Qout=trapz(mesh.t,medium.Iconst.*abs(EtF).^2).*(medium.area_hcf/2);
            Qout2=trapz(mesh.t,medium.Iconst.*abs(EtS).^2).*(medium.area_hcf/2);
            Qerror=(Qout-Qout2)/Qout;
%                         myplot1('intensity','no','no',mesh.f,EfF)
%                         xlim([-0e15 1e15])
        %                 myplot1('field','envelope','no',mesh.t,EtF.*exp(pulse1.carrier),EtS.*exp(pulse2.carrier))    
        %                 xlim([-100e-15 200e-15])
            plot(mesh.t,[real(EtF);abs(EtF)])
            xlim([-100e-15 100e-15])

            [~,x_SHInitmid]=findpeaks(abs(pulse2.Ert),mesh.t,'NPeaks',1,'SortStr','descend');
            [~,x_SHmid]=findpeaks(abs(EtS),mesh.t,'NPeaks',1,'SortStr','descend');
            dt=abs(x_SHInitmid-x_SHmid);
            title(['z=',num2str(m*mesh.dz),'  ','Q1=',num2str(Qout),'  ','Q2=',num2str(Qout2),'  ','Qin=',num2str(pulse1.Energy+pulse2.Energy),'  ','dt=',num2str(dt)]);                                  
            pause(0.1)   
        %% Ionization Level
        [n_e,Eg,n_gas]=calc_eDensityADK((EtF),mesh,medium,pulse1,beam);
        IonizLvl=[IonizLvl,max(n_e)/n_gas];
    end
    
    
end
    
function [newErf]=do_oneRKstep(mesh,pulse1,pulse2,beam,medium,Ef1,Ef2,m,h)
%calculate next step via Runge Kutta
       k1= mesh.dz.*calcfunction(mesh,pulse1,pulse2,beam,medium,Ef1,Ef2,1);   
       k2 = mesh.dz.*calcfunction(mesh,pulse1,pulse2,beam,medium,Ef1+k1./2,Ef2,2);        
       k3 = mesh.dz.*calcfunction(mesh,pulse1,pulse2,beam,medium,Ef1+k2./2,Ef2,2);           
       k4 = mesh.dz.*calcfunction(mesh,pulse1,pulse2,beam,medium,Ef1+k3,Ef2,2);
       newErf=Ef1 + (k1+2.*k2+2.*k3+k4)./6;   
% plot(mesh.t,[abs(k1);abs(k2);abs(k3);abs(k4)])       
end

function [result]=calcfunction(mesh,pulse1,pulse2,beam,medium,Ef1,Ef2,m)
%calculate diff. equation for Runge Kutta
        Et1=myifft(Ef1,mesh);
        Et2=myifft(Ef2,mesh);
    %SPM        
        It1=medium.Iconst.*abs(Et1).^2;
        SPM=-1i*(medium.n2*(pulse1.w0)./const.c).*It1;
    %XPM
        It2=medium.Iconst.*abs(Et2).^2;
        XPM=-1i*(2*medium.n2*pulse1.w0/const.c).*It2;%refractive index is 2 times higher for XPM!
    %Crossterm
%         CRS=-1i*(medium.n2*(pulse1.w0)./const.c).*medium.Iconst.*(Et1.*conj(Et2));
    %SST
        SPMSST2=(2*pi.*(mesh.f))./(pulse1.w0).*myfft(SPM.*Et1,mesh);%
%         SPMSST2(1:mesh.fmid)=0;
        LRbounds=find_bounds(SPMSST2);
        fwidth=mesh.df.*(LRbounds(1,2)-LRbounds(1,1));
%         smoothfct=calc_supergaussian(mesh.f,fwidth*2,10,0);%pulse1.w0/(2*pi)
        smoothfct=calc_smartfilter(mesh,SPMSST2,pulse1);
    %Ionization
        [n_e,Eg,n_gas]=calc_eDensityADK((Et1),mesh,medium,pulse1,beam);
    %Ionization Loss
%         dNedt=zeros(1,mesh.flength); dNedt(1:(end-1))=diff(n_e)/mesh.dt;
        dNedt=gradient(n_e,mesh.dt);%diff(n_e)/mesh.dt;  (1:(end-1))
        smooth=calc_supergaussian(mesh.t,beam.delT.*2,10,0);
%         It=medium.Iconst.*abs(Et1).^2;
        It=medium.Iconst.*(real(Et1)).^2;
        ION=-Eg.*(dNedt)./(2.*It);
        ION(isnan(ION))=0;  
        ION(isinf(ION))=0;
%         winst=-[diff(abs(SPM))./mesh.dt,0];    %instantaneous frequency
    %Ionization BlueShift
        wp2=const.e^2.*n_e./(const.eps0*const.m_e);
        wp2=wp2.*sqrt(3);
        [dn_p,nL]=calc_RefrIndexPlasma(mesh,beam,medium,n_e,n_gas);%refr. Index of Plasma if Plasma+Gas mixed
        PLSM=(1i*pulse1.w0/const.c).*wp2./(2.*pulse1.w0.^2);
        PLSM2=(1i*pulse1.w0/const.c).*dn_p;
        PLSM3=(1i*pulse1.w0/const.c).*(1-sqrt(1-wp2./(pulse1.w0.^2)));

result=myfft((0.*XPM+0*SPM+0.*ION+0.*PLSM3).*Et1,mesh)+1*SPMSST2.*smoothfct;%myfft((SPM).*Et1,mesh);
result(1:mesh.fmid)=0;
end