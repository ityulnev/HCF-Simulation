%Create structure with all needed parameters specific to the beam
classdef beam_init

    properties
    wavelength,f0,w0,Q_In,Q_In2,Q_Out,delT,tau,tau2,peakpower,n_cycles,Ipeak,tfwhm
    end
    
    methods 
        function s=beam_init
        s.wavelength=800e-9;%[m]
        s.f0=const.c/s.wavelength;
        s.w0=2*pi*s.f0;
        s.n_cycles=1;
        s.Q_In=2.2e-3.*(s.n_cycles/1);%[J]
        s.Q_In2=0.1e-3;%[J]
        s.Q_Out=2.2e-3;%[J]
        s.Ipeak=1;%2*1e19;%[W/m^2] Peak Intensity 
        %% Pulse Duration
        s.delT=35e-15;
%         s.delT=s.n_cycles/s.f0;%[s] Pulse duration
        end
    end
end