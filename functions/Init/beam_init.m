%Create structure with all needed parameters specific to the beam
classdef beam_init

    properties
    wavelength,f0,w0,Q_In,Q_In2,Q_Out,t_fwhm,tau,tau2,peakpower
    end
    
    methods 
        function s=beam_init
        s.wavelength=800e-9;%[m]
        s.f0=const.c/s.wavelength;
        s.w0=2*pi*s.f0;
        s.Q_In=2.0e-3;%[J]
        s.Q_In2=0.1e-3;%[J]
        s.Q_Out=1e-3;%[J]
        s.t_fwhm=35e-15;%[s]
%         s.t_fwhm2=s.t_fwhm/sqrt(s.order);
%         s.tau=s.t_fwhm/(2*sqrt(log(2)));
%         s.tau2=s.t_fwhm2/(2*sqrt(log(2)));
%         s.peakpower=s.Q_In/s.t_fwhm;%[J/s]
        end
    end
end