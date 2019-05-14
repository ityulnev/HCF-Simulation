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
        s.Q_In=1.0e-3;%[J]
        s.Q_In2=1e-3;%[J]
        s.Q_Out=1e-3;%[J]
        s.t_fwhm=7e-15;%[s]
        end
    end
end