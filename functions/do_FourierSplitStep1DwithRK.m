%Calculate propagated electric field after distance Lz via fourier split
%step
%d/dz(E) from SVEA in 1D
%Implemented Effects: GVD, SPM, Attenuation, Ionization
function [Et,Ef]=do_FourierSplitStep1DwithRK(mesh,beam,fiber,pulse)
h=mesh.dz/2;

% Group Velocity
% dndw=[diff(fiber.n)./(mesh.df*2*pi),0];
% dndw_ext=[zeros(1,mesh.fbound-1),dndw];
% 
% kph=fiber.n_ext./const.c;
% kgr=kph+dndw_ext.*((2*pi.*mesh.f)./const.c);

         

Et=pulse.Ert;

    for m=1:(mesh.zlength)
            
            Et=do_oneRKstep(mesh,pulse,beam,fiber,Et,m);
            Ef=myfft(Et,mesh.dt);
            
%             %K
%             KK=1i.*fiber.k;         
%             %k0
%             kk=1i.*fiber.k0;             
%             %group Velocity GV
%             GVV=(1i.*((mesh.w-beam.w0).*(kgr-kph)));
%             GV=1i.*(2*pi.*mesh.f-beam.w0).*kgr;
%             PV=1i.*(mesh.w-beam.w0).*kph;
            
            
            %Group velocity dispersion
            GVD=1i.*((2*pi.*mesh.f)+beam.w0).^2.*fiber.Beta2/2;
%             GVD(1:mesh.fbound)=0;
                       
            Ef=Ef.*exp((GVD).*h).*exp(-fiber.alpha*h);                      
            Et=myifft(Ef,mesh.dt);
            
            Qin=pulse.Energy;
            Qout=trapz(mesh.t,fiber.Iconst.*abs(Et).^2).*fiber.area_hcf;
%             plot(mesh.f,[abs(myifft(pulse.Ert,mesh.dt)).^2;abs(Ef).^2])
%             xlim([2.5e14 5e14])
%             plot(mesh.t,[abs(pulse.Ert).^2;abs(Et).^2])
            plot(mesh.t,[real(pulse.Ert);abs(pulse.Ert);real(Et);abs(Et)]);
            xlim([-60e-15 60e-15])

title(['z=',num2str(m*mesh.dz),'  ','Qout=',num2str(Qout)]);                                  
pause(0.1)                      
    end
    
end
    
    function [newErt]=do_oneRKstep(mesh,pulse,beam,fiber,Ert,m)
%calculate next step via Runge Kutta
h=mesh.dz/2;
       k1= h.*calcfunction(mesh,pulse,beam,fiber,Ert);   
       k2 = h.*calcfunction(mesh,pulse,beam,fiber,Ert+k1./2);        
       k3 = h.*calcfunction(mesh,pulse,beam,fiber,Ert+k2./2);           
       k4 = h.*calcfunction(mesh,pulse,beam,fiber,Ert+k3);
       newErt=Ert + (k1+2.*k2+2.*k3+k4)./6;       
    end

function [result]=calcfunction(mesh,pulse,beam,fiber,Et)
%calculate diff. equation for Runge Kutta
    %SPM        
        It=fiber.Iconst.*abs(Et).^2;
        SPM=1i*(fiber.n2*beam.w0/const.c).*It;
    %SST
        SPMSST2=(2*pi.*(mesh.f-beam.f0))./(beam.w0).*myfft(SPM.*Et,mesh.dt);%
        SPMSST2(abs(SPMSST2)<max(abs(SPMSST2))/(1e6))=0;%take care of numerical errors which might occur
        SPMSST=myifft(SPMSST2,mesh.dt);
    %Ionization
        [n_e,Eg]=calc_eDensityADK(abs(Et),mesh,fiber);
        dNedt=zeros(1,mesh.flength);
        dNedt(1:(end-1))=diff(n_e)/mesh.dt;          
        ION=-Eg.*(dNedt)./(2.*It);
        ION(isnan(ION))=0;                        
result=(0.*ION).*Et+SPMSST;
end