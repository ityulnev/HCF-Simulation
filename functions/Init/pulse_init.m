%Create structure that contains properties of the pulse in t,f and r domain
classdef pulse_init
    
   properties
   w0,t_fwhm,tau0,t0,ef,et,er,Ert,carrier,A0,A0f,Fluence,Energy,Erf,Energyf,order,Irf,Irt,ptmid,pfmid
   end
methods
    function s=pulse_init(mesh,beam,medium,Q_In,t0,order)    
        s.order=order;
        s.t_fwhm=beam.t_fwhm/sqrt(s.order);% Pulse Duration @ FWHM
        s.tau0=s.t_fwhm/(2*sqrt(log(2)));%tau0 for Gaussian variance
        s.w0=beam.w0.*s.order;
        s.t0=t0;%Time Delay in s
        timedelay=(s.t0*1i*2*pi.*(mesh.f));     
        s.carrier=1i*s.w0.*mesh.t;% Oscillation of carrier wave
        s.ef=exp(-(2*pi.*(mesh.f)).^2.*s.tau0^2./2-timedelay);
        et=myifft(s.ef,mesh);
        s.et=et.*exp(1i.*s.w0.*(-s.t0));%.*exp(s.carrier); NO INITIAL CARRIER ### But shift of carrier from time delay
        s.Fluence=Q_In/medium.area_hcf;
        s.A0=sqrt(s.Fluence/(sum(medium.Iconst.*abs(s.et).^2)*mesh.dt)); %
        s.A0f=sqrt(s.Fluence/(sum(medium.Iconst.*abs(s.ef).^2)*mesh.df)); 
        n_gaussian=1;
        s.er=exp(-(((mesh.r).^2./((medium.r_mode)^2))).^n_gaussian);
%% Field and Intensity
        s.Ert=s.A0.*transpose(s.er).*s.et;
%         s.Erf=s.A0f.*transpose(s.er).*s.ef; 
        s.Erf=myfft(s.Ert,mesh);
        s.Irt=medium.Iconst.*abs(s.Ert).^2;
        s.Irf=medium.Iconst.*abs(s.Erf).^2;
        
        s.ptmid=find(max(s.Irt)==s.Irt);
        s.pfmid=find(max(s.Irf)==s.Irf);       
%% Calculate Pulse energy        
        if ndims(s.Ert)>2
            s.Energy=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.t,s.Irt,2),1);
            s.Energyf=2*pi.*trapz(mesh.r,transpose(mesh.r).*trapz(mesh.f,s.Irf,2),1);
        else
            s.Energy=medium.area_hcf.*trapz(mesh.t,s.Irt); 
            s.Energyf=medium.area_hcf.*trapz(mesh.f,s.Irf); 
        end             
%% Test: Energy conservation & FWHM check for FFT
        tolerance=1e-4;
        if abs(s.Energyf-s.Energy)/s.Energy <tolerance
            %do nothing  
        else
            warning('pulse_init: Energy of Et and Ef not conserved!')    
        end
        
        fwhmF=calc_fwhm(mesh.f,s.Irf);
        fwhmT=calc_fwhm(mesh.t,s.Irt);
        tolerance=2e-2;
        if abs(fwhmT*fwhmF-0.44)/0.44 <tolerance
            %do nothing  
        else
            warning('pulse_init: FWHM of Et and Ef not conserved!')    
        end
                
    end
end
end