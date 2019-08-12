%Calculate propagated electric field after distance Lz via fourier split
%step
%d/dz(E) from SVEA in 1D
%###################VALID ONLY FOR PHASE MATCHED E1 and E2##########
function [EtF]=do_2Color1DwithSummedFields(mesh,beam,fiber,pulse1,pulse2)
h=mesh.dz;
EtF=pulse1.Ert+pulse2.Ert;
    for m=1:(mesh.zlength)
            %% Self phase modulation SPM and XPM
            newEtF=do_oneRKstep(mesh,pulse1,pulse2,beam,fiber,EtF,pulse1.Ert,m,h);
            EfF=myfft(newEtF,mesh);                
            %% Group velocity dispersion GVD
            %Fundamental
            GVD=-1i.*(2*pi.*mesh.f).^2.*fiber.kGVD(pulse1.pfmid)./2; 
%             GVD=0;
            EfF=EfF.*exp((GVD).*h).*exp(-fiber.alpha/2*h);    
            EfF(isinf(EfF))=0;
            EtF=myifft(EfF,mesh);
            %% Plot
    Qout=trapz(mesh.t,fiber.Iconst.*abs(EtF).^2).*fiber.area_hcf;
%                 myplot1('intensity','no','no',mesh.f,EfF)
%                 xlim([-1e15 1e15])
                myplot1('field','envelope','no',mesh.t,EtF.*exp(pulse1.carrier))    
                xlim([-100e-15 200e-15])                 
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
        SPMSST2=(2*pi.*(mesh.f)+pulse1.w0)./(pulse1.w0).*myfft((SPM).*Et1,mesh);%
        LRbounds=find_bounds(SPMSST2);
        fwidth=mesh.df.*(LRbounds(1,2)-LRbounds(1,1));
        smoothfct=calc_supergaussian(mesh.f,fwidth.*3,10);
        SPMSST=myifft(SPMSST2.*smoothfct,mesh);                     
   %Ionization
%         [n_e,Eg]=calc_eDensityADK(abs(Et1),mesh,fiber);
%         dNedt=zeros(1,mesh.flength);
%         dNedt(1:(end-1))=diff(n_e)/mesh.dt;
%         It=fiber.Iconst.*abs(Et1).^2;
%         ION=-Eg.*(dNedt)./(2.*It);
%         ION(isnan(ION))=0;  
%         winst=[diff(imag(SPM))./mesh.dt,0];    %instantaneous frequency
    
result=(SPM.*0).*Et1+SPMSST;
end