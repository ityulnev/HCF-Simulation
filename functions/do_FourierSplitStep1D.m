%Calculate propagated electric field after distance Lz via fourier split
%step
%d/dz(E) from SVEA in 1D
%Implemented Effects: GVD, SPM, Attenuation, Ionization
function [Et,Ef]=do_FourierSplitStep1D(mesh,beam,fiber,pulse)
h=mesh.dz/2;

% Group Velocity
% dndw=[diff(fiber.n)./(mesh.df*2*pi),0];
% dndw_ext=[zeros(1,mesh.fbound-1),dndw];
% 
% kph=fiber.n_ext./const.c;
% kgr=kph+dndw_ext.*((2*pi.*mesh.f)./const.c);

         

Et=pulse.Ert;

    for m=1:(mesh.zlength)
            
            It=fiber.Iconst.*abs(Et).^2;
            %SPM
            SPM=1i*(fiber.n2*beam.w0/const.c).*It;
            %SST
            SPMSST2=(2*pi.*(mesh.f+beam.f0))./(beam.w0).*myifft(SPM,mesh.dt);
            SPMSST2(abs(SPMSST2)<1e-19)=0;%take care of numerical errors
            SPMSST=myfft(SPMSST2,mesh.dt);
            %Ionization
            [n_e,Eg]=calc_eDensityADK(abs(Et),mesh,fiber);
            dNedt=zeros(1,mesh.flength);
            dNedt(1:(end-1))=diff(n_e)/mesh.dt;          
            ION=-Eg.*(dNedt)./(2.*It);
            ION(isnan(ION))=0;

            Et=Et.*exp((SPMSST+ION).*h).*exp(-fiber.alpha*h);
            Ef=myifft(Et,mesh.dt);
            
%             %K
%             KK=1i.*fiber.k;         
%             %k0
%             kk=1i.*fiber.k0;             
%             %group Velocity GV
%             GVV=(1i.*((mesh.w-beam.w0).*(kgr-kph)));
%             GV=1i.*(2*pi.*mesh.f-beam.w0).*kgr;
%             PV=1i.*(mesh.w-beam.w0).*kph;
                        
            %Group velocity dispersion
            GVD=1i.*((2*pi.*mesh.f)-beam.w0).^2.*fiber.Beta2/2;
            GVD(1:mesh.fbound)=0;            
            
            Ef=Ef.*exp((GVD).*h);                      
            Et=myfft(Ef,mesh.dt);
            
            Qin=pulse.Energy;
            Qout=trapz(mesh.t,fiber.Iconst.*abs(Et).^2).*fiber.area_hcf;
            plot(mesh.f,[abs(myifft(pulse.Ert,mesh.dt)).^2;abs(Ef).^2])
            xlim([2.5e14 5e14])
%             plot(mesh.t,[abs(pulse.Ert).^2;abs(Et).^2])
%             plot(mesh.t,[real(pulse.Ert);abs(pulse.Ert);real(Et);abs(Et)]);
%             xlim([-50e-15 50e-15])

title(['z=',num2str(m*mesh.dz),'  ','Qout=',num2str(Qout)]);                                  
pause(0.1)   

                   
    end
end
