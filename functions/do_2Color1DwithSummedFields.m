%Calculate electric field propagation with steps in 0:dz:Lz
%Electric field assumed to be without carrier and sum of 2 pulses
%d/dz(E) from SVEA in 1D
%Implemented Effects: GVD, SPM+SST, XPM, Attenuation, Ionization
%###################VALID ONLY FOR PHASE MATCHED E1 and E2##########
function [EtF,EtF1,EtS1]=do_2Color1DwithSummedFields(mesh,beam,medium,pulse1,pulse2)
h=mesh.dz;
%% Initial chirp
It1=medium.Iconst.*abs(pulse1.Ert).^2;
SPM1=-1i*(medium.n2*(pulse1.w0)./const.c).*It1;
EtF1=(pulse1.Ert).*exp(SPM1.*0);
%GVD chirp
EtF1=myifft(pulse1.Erf.*exp(-1i*20.*(2*pi.*mesh.f-pulse1.w0).^2.*medium.kGVD(pulse1.pfmid)./2),mesh);
EtF1=pulse1.Ert;
Peakintensityscale=sqrt(0.01);
EtS1=myifft(myfft(EtF1,mesh).*exp(-1i.*(2*pi.*(mesh.f)).*pulse2.t0),mesh).*Peakintensityscale;
Etin=(EtF1+EtS1);
EtF=Etin;

%% No chirp
% EtF=(pulse1.Ert)+(pulse2.Ert);
    for m=1:(mesh.zlength)              
            %% Self phase modulation SPM and XPM
            newEtF=do_oneRKstep(mesh,pulse1,pulse2,beam,medium,EtF,pulse1.Ert,m,h);
            EfF=myfft(newEtF,mesh);  
        %% Group velocity dispersion GVD
            GVD=-1i.*(2*pi.*mesh.f-pulse1.w0).^2.*medium.kGVD(pulse1.pfmid)./2; 
%             GVD=0;
            EfF=EfF.*exp((GVD).*h).*exp(-beam.alpha/2*h);    
            EfF(isinf(EfF))=0;
            EfF(1:mesh.fmid)=0;
            EtF=myifft(EfF,mesh);
            %% Plot
    Qout=trapz(mesh.f,medium.Iconst.*abs(EfF).^2).*medium.area_hcf;
%                 myplot1('intensity','no','no',mesh.f,myfft(EtF,mesh),myfft(pulse1.Ert+pulse2.Ert,mesh))
%                 xlim([-1e15 1e15])
                myplot1('field','envelope','no',mesh.t,EtF)    
                xlim([-200e-15 200e-15])                 
    title(['z=',num2str(m*mesh.dz),'  ','Qout=',num2str(Qout),'  ','Qin=',num2str(pulse1.Energy+pulse2.Energy)]);                                  
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
%         It2=fiber.Iconst.*abs(Et2).^2;
%         XPM=-1i*(2*fiber.n2*pulse1.w0/const.c).*It2;%refractive index is 2 times higher for XPM!
    %SST
%         SPMSST2=(2*pi.*(mesh.f)+pulse1.w0)./(pulse1.w0).*myfft((SPM).*Et1,mesh);%
%         LRbounds=find_bounds(SPMSST2);
%         fwidth=mesh.df.*(LRbounds(1,2)-LRbounds(1,1));
%         smoothfct=calc_supergaussian(mesh.f,fwidth.*3,10);
%         SPMSST=myifft(SPMSST2.*smoothfct,mesh);                     
   %Ionization
%         [n_e,Eg]=calc_eDensityADK(abs(Et1),mesh,fiber);
%         dNedt=zeros(1,mesh.flength);
%         dNedt(1:(end-1))=diff(n_e)/mesh.dt;
%         It=fiber.Iconst.*abs(Et1).^2;
%         ION=-Eg.*(dNedt)./(2.*It);
%         ION(isnan(ION))=0;  
%         winst=[diff(imag(SPM))./mesh.dt,0];    %instantaneous frequency
    
result=(1.*SPM).*Et1;%+1.*SPMSST;
end