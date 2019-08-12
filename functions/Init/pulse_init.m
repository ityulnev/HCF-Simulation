%Create structure that contains properties of the pulse in t,f and r domain
classdef pulse_init
    
   properties
   w0;
   t_1sigma,r_1sigma,t0,ef,et,er,Ert,carrier,A0,A0f,Fluence,Energy,Erf,Energyf,order,Irf,Irt,ptmid,pfmid,CEP,cep,fwhmF,fwhmT,tau02,Ipeak,IpeakTheo,PpeakTheo
   end
methods
    function s=pulse_init(mesh,beam,medium,Q_In,t0,order,cep)    
        %% Set correct pulse duration and calculate carrier and timedelay phase
        s.order=order;
        % t_1sigma = dT_FW/2 @Imax/2
        s.t_1sigma=beam.delT/(2*sqrt(log(2)));                                
        % t_1sigma = dT_FW/2 @Imax/e^2
%         s.t_1sigma=beam.delT/(2*sqrt(2));  
        
        s.t_1sigma=s.t_1sigma/sqrt(s.order);
        s.w0=beam.w0.*s.order;
        s.t0=t0;%Time Delay in s
        tdelayphase=1i.*(2*pi.*(mesh.f)).*s.t0;     
        s.carrier=1i*s.w0.*mesh.t;% Oscillation of carrier wave
        efold=exp(-(2*pi.*(mesh.f)).^2.*s.t_1sigma^2./2-tdelayphase);
        efmirror=-exp(-(2*pi.*(mesh.f+beam.f0*2)).^2.*s.t_1sigma^2./2-tdelayphase);
        %% This filter takes care of negative frequencies for very broadband pulses approaching f=0
        %Case: 0.5*F_fwhm>[0,f0]
        tanhfilter=calc_tanhfilter(2*pi.*mesh.t.*(mesh.f+beam.f0),mesh.f0equal0);
%         tanhfilter=-tanh(2*pi*mesh.t.*(mesh.f+beam.f0));
%         tanhmax=find(max(tanhfilter)==tanhfilter,1);
%         tanhfilter(tanhmax:end)=1;
%         tanhfilter(1:mesh.fequal0)=0;
        s.ef=efold.*tanhfilter;
        %% Gaussian mirror filter
%         s.ef=efold+efmirror;
%         s.ef=efold;
%         s.ef(1:mesh.fequal0)=0;
        %%
        et=myifft(s.ef,mesh);
        s.cep=cep;
        s.CEP=1i.*s.w0.*(-s.cep*pi./s.w0);% Carier envelope phase phi=w0*t1
        s.et=et.*exp(1i.*s.w0.*(-s.t0)).*exp(s.CEP).*exp(s.carrier);% NO INITIAL CARRIER ### But shift of carrier from time delay
        if beam.Ipeak==1 
            s.Fluence=Q_In/(medium.area_hcf/2);
            s.A0=sqrt(s.Fluence/(sum(medium.Iconst.*abs(s.et).^2)*mesh.dt)); %
%             s.A0f=sqrt(s.Fluence/(sum(medium.Iconst.*abs(s.ef).^2)*mesh.df)); 
        else
            s.A0=sqrt(beam.Ipeak/medium.Iconst)./max(abs(s.et));
        end
        n_gaussian=1;
        s.r_1sigma=medium.r_mode/sqrt(2);%1sigma variance of gaussian from beam radius @I/e^2  
        s.er=exp(-(((mesh.r).^2./(2*(s.r_1sigma)^2))).^n_gaussian);
%% Field and Intensity
        s.Ert=s.A0.*transpose(s.er).*s.et;
%         s.Erf=s.A0f.*transpose(s.er).*s.ef; 
        s.Erf=myfft(s.Ert,mesh);
        s.Irt=medium.Iconst.*abs(s.Ert).^2;
        s.Irf=medium.Iconst.*abs(s.Erf).^2;
        
        s.ptmid=find(max(s.Irt)==s.Irt);
        s.pfmid=find(max(s.Irf)==s.Irf);       
%%      Peak Powe&Intensity
        s.Ipeak=max(s.Irt);
        s.PpeakTheo=0.94*(Q_In/(beam.delT/sqrt(order)));
        s.IpeakTheo=s.PpeakTheo/(medium.area_hcf/2);%Peak intensity for Gaussian beam from peak power
        %*sqrt(log(2)/2)
        %Test
        if beam.Ipeak==1
        tolerance=1e-2;
            if abs(s.IpeakTheo-s.Ipeak)/s.IpeakTheo <tolerance
                %do nothing  
            else
                warning('pulse_init: Peak Intensity of Irt deviates from theory!')    
            end
        end
        %% Calculate Pulse energy        
        if ndims(s.Ert)>2
            s.Energy=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t,s.Irt,2),1);
            s.Energyf=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.f,s.Irf,2),1);
        else
            s.Energy=medium.area_hcf.*trapz(mesh.t,s.Irt); 
            s.Energyf=medium.area_hcf.*trapz(mesh.f,s.Irf); 
        end             
%% Test: Energy conservation & FWHM check with FFT
        tolerance=1e-4;
        if abs(s.Energyf-s.Energy)/s.Energy <tolerance
            %do nothing  
        else
            warning('pulse_init: Energy of Et and Ef not conserved!')    
        end
        s.fwhmF=calc_fwhm(mesh.f,s.Irf);
        s.fwhmT=calc_fwhm(mesh.t,s.Irt);
        tolerance=2e-2;
        tbp=2*log(2)/pi;
        if abs(s.fwhmT*s.fwhmF-tbp)/tbp <tolerance
            %do nothing  
        else
            warning('pulse_init: FWHM of Et and Ef not conserved!')    
        end
                
    end
end
end