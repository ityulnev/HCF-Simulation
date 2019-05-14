%Calculate propagated electric field after distance Lz via fourier split
%step
%d/dz(E) from SVEA in 1D
%Implemented Effects: GVD, SPM+SST, XPM, Attenuation, Ionization
function [EtF,EtS]=do_2ColorFourierSplitStep1DwithRK(mesh,beam,fiber,pulse1,pulse2)
h=mesh.dz;%Half step size for fourier split step

index2f0=find(mesh.f==(pulse2.w0/(2*pi)-beam.f0),1);
% Lwalkoff=pulse2.t_fwhm/abs(fiber.kGV(indexf0)-fiber.kGV(index2f0)); % Walkoff length for 2 pulses copropagating
Twalkoff=abs(fiber.kGV(pulse1.pfmid)-fiber.kGV(index2f0)).*mesh.Lz;
EtF=pulse1.Ert;% Fundamental pulse
EtS=pulse2.Ert;% Second Harmonic pulse
    for m=1:(mesh.zlength)
            %% Self phase modulation SPM and XPM
            newEtF=do_oneRKstep(mesh,pulse1,pulse2,beam,fiber,EtF,EtS,m,h);
            EfF=myfft(newEtF,mesh); 
            newEtS=do_oneRKstep(mesh,pulse2,pulse1,beam,fiber,EtS,EtF,m,h);
            EfS=myfft(newEtS,mesh);           
%             EfF=myfft(EtF,mesh); 
%             EfS=myfft(EtS,mesh);               
            %% Group Velocity difference        
%             GV1=-1i.*(2*pi.*mesh.f-pulse1.w0).*fiber.kGV(indexf0);
%             GV2=-1i.*(2*pi*mesh.f-pulse2.w0).*fiber.kGV(index2f0);
%             dGV=-1i.*(2*pi*mesh.f-pulse2.w0).*(fiber.kGV(index2f0)-fiber.kGV(indexf0));
            dGV=-1i.*(2*pi*mesh.f).*(fiber.kGV(index2f0)-fiber.kGV(pulse1.pfmid));

            %% Group velocity dispersion GVD
            %Fundamental
            GVD=-1i.*(2*pi.*mesh.f).^2.*fiber.kGVD(pulse1.pfmid)./2; 
%             GVD=0;
%             TOD=-1i.*((2*pi.*mesh.f)-pulse1.w0).^3.*fiber.kTOD(indexf0)./6; 

            EfF=EfF.*exp((GVD).*h).*exp(-fiber.alpha/2*h);    
            EfF(isinf(EfF))=0;
            EtF=myifft(EfF,mesh);
            %Second Harmonic
            GVD2=-1i.*(2*pi.*mesh.f).^2.*fiber.kGVD(index2f0)./2;
%             GVD2=0;
%             TOD2=-1i.*((2*pi.*mesh.f)-pulse2.w0).^3.*fiber.kTOD(index2f0)./6;
            EfS=EfS.*exp((GVD2+dGV).*h).*exp(-fiber.alpha/2*h);% -> Energy loss                     
            EfS(isinf(EfS))=0;
            EtS=myifft(EfS,mesh);
            %% Plot
    Qout=trapz(mesh.t,fiber.Iconst.*abs(EtF).^2).*fiber.area_hcf;
    Qout2=trapz(mesh.t,fiber.Iconst.*abs(EtS).^2).*fiber.area_hcf;
    Qerror=(Qout-Qout2)/Qout;
%                 myplot1('intensity','no','no',mesh.f,EfF,EfS)
%                 xlim([-1e15 1e15])
                myplot1('field','envelope','no',mesh.t,EtF.*exp(pulse1.carrier),EtS.*exp(pulse2.carrier))    
                xlim([-100e-15 200e-15])

    [~,x_SHInitmid]=findpeaks(abs(pulse2.Ert),mesh.t,'NPeaks',1,'SortStr','descend');
    [~,x_SHmid]=findpeaks(abs(EtS),mesh.t,'NPeaks',1,'SortStr','descend');
    dt=abs(x_SHInitmid-x_SHmid);
    title(['z=',num2str(m*mesh.dz),'  ','Q1=',num2str(Qout+Qout2),'  ','Qin=',num2str(pulse1.Energy+pulse2.Energy),'  ','dt=',num2str(dt)]);                                  
    pause(0.1)                      
    end
    
end
    
    function [newErt]=do_oneRKstep(mesh,pulse1,pulse2,beam,fiber,Et1,Et2,m,h)
%calculate next step via Runge Kutta
       k1= mesh.dz.*calcfunction(mesh,pulse1,pulse2,beam,fiber,Et1,Et2,1);   
       k2 = mesh.dz.*calcfunction(mesh,pulse1,pulse2,beam,fiber,Et1+k1./2,Et2,2);        
       k3 = mesh.dz.*calcfunction(mesh,pulse1,pulse2,beam,fiber,Et1+k2./2,Et2,2);           
       k4 = mesh.dz.*calcfunction(mesh,pulse1,pulse2,beam,fiber,Et1+k3,Et2,2);
       newErt=Et1 + (k1+2.*k2+2.*k3+k4)./6;   
% plot(mesh.t,[abs(k1);abs(k2);abs(k3);abs(k4)])       
    end

function [result]=calcfunction(mesh,pulse1,pulse2,beam,fiber,Et1,Et2,m)
%calculate diff. equation for Runge Kutta
        
    %SPM        
        It1=fiber.Iconst.*abs(Et1).^2;
        SPM=-1i*(fiber.n2*(pulse1.w0)./const.c).*It1;
    %XPM
        It2=fiber.Iconst.*abs(Et2).^2;
        XPM=-1i*(2*fiber.n2*pulse1.w0/const.c).*It2;%refractive index is 2 times higher for XPM!
    %Crossterm
        CRS=-1i*(fiber.n2*(pulse1.w0)./const.c).*fiber.Iconst.*(Et1.*conj(Et2));
    %SST
%         SPMSST2=(2*pi.*(mesh.f)+pulse1.w0)./(pulse1.w0).*myfft((SPM).*Et1,mesh);%
%         LRbounds=find_bounds(SPMSST2);
%         fwidth=mesh.df.*(LRbounds(1,2)-LRbounds(1,1));
%         smoothfct=calc_supergaussian(mesh.f,fwidth.*3,10);
%         SPMSST=myifft(SPMSST2.*smoothfct,mesh);                     
%    %Ionization
%         [n_e,Eg]=calc_eDensityADK(abs(Et1),mesh,fiber);
%         dNedt=zeros(1,mesh.flength);
%         dNedt(1:(end-1))=diff(n_e)/mesh.dt;
%         It=fiber.Iconst.*abs(Et1).^2;
%         ION=-Eg.*(dNedt)./(2.*It);
%         ION(isnan(ION))=0;  
%         winst=[diff(imag(SPM))./mesh.dt,0];    %instantaneous frequency
    
result=(XPM+SPM+CRS).*Et1;
end