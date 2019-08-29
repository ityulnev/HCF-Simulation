%Create structure with all needed parameters specific to the beam
classdef beam_init

    properties
    wavelength,f0,w0,Q_In,Q_In2,Q_Out,delT,tau,tau2,peakpower,n_cycles,Ipeak,tfwhm,alpha
    end
    
    methods 
        function s=beam_init
        s.wavelength=800e-9;%[m]
        s.f0=const.c/s.wavelength;
        s.w0=2*pi*s.f0;
        s.n_cycles=3;
        s.Q_In=2e-3.*(s.n_cycles/3);%[J]
        s.Q_In2=0.1e-3;%[J]
        s.Q_Out=2e-3;%[J]
        s.alpha=log(s.Q_In/s.Q_Out);% 0.8 attenuation coeffiecient
        s.alpha=0;
     
        s.Ipeak=1;%[W/m^2] Peak Intensity 
        %% Pulse Duration
        s.delT=35e-15;
%         s.delT=s.n_cycles/s.f0;%[s] Pulse duration
        end
    end
end