%Calculates W_adk for a given E field and a given Gas
function [n_e,Eg,n_gas]=calc_eDensityADK(Et,mesh,fiber,pulse,beam)
% Et=0.5.*(Et+conj(Et));
Et=real(Et);
smooth=calc_supergaussian(mesh.t,beam.delT.*2,10,0);
Et=smooth.*abs(Et);

%          plot(mesh.t,Et)
%          pause(0.1)
%Quantum Numbers
l=1;
m=0;
%Gastype
switch fiber.gas
    case 'Neon'
            l=1;
            m=0;
            Eg=21.565*const.e;                                             %[J] from http://www.periodensystem.info/elemente/neon/
            n_gas=2.686e25*fiber.pressure;                                                %1/m^3 
    case 'Argon'
            Eg=15.76*const.e;                                              %[J] from http://www.periodensystem.info
            n_gas=2.7e25*fiber.pressure;                                                  %1/m^3   
    case 'Xenon'
            Eg=12.13*const.e;                                              %[J] from http://www.periodensystem.info
            n_gas=2.4e25*fiber.pressure;                                                  %1/m^3  
    case 'Helium'
            l=1;
            m=0;
            Eg=24.587*const.e;
            n_gas=2.6856e25*fiber.pressure;
            
end
%%
f=(2*l+1)*factorial(l+abs(m))/(2^abs(m)*factorial(abs(m))*factorial(l-m));
Z=1;
%% Conversion into Atomic Units
rbohr=4*pi*const.eps0*const.hbar^2/(const.m_e*const.e^2);   
E_hartree=const.m_e*const.e^4/(4*const.eps0^2*const.h^2);                  %Hartree Energy [J]
Efield_hartree=E_hartree/(rbohr*const.e);                                  %Efield from Hartree Energy [V/m]
Eg_au=Eg/E_hartree;                                                        %Energy into atomic units a.u.
E_laser_au=Et./Efield_hartree;                                        %Field into atomic untis a.u.  
dt_atu=mesh.dt./(const.hbar/E_hartree);                             %Time into atomic time units a.t.u.

%% Calc ionization density n_e with ADK model
n_star=Z*(1/sqrt(2*Eg_au));                                                %Effective principal quantum number
E_0_au=(2*Eg_au)^(3/2);
fraction=E_0_au./E_laser_au;
exponent=exp(-(2/3).*fraction);
w_factor=(E_laser_au).^(1.5-2*n_star);
w_factor(isinf(w_factor))=0;
prod=w_factor.*exponent;
prod(isinf(prod))=0;
prod(isnan(prod))=0;
W_adk=f/(8*pi*n_star).*((4*exp(1)*E_0_au./(n_star)).^(2*n_star)).*sqrt(3./(pi*E_0_au*(2*Eg_au))).*prod;  %Ionization rate
P_adk=1-exp(cumsum(-dt_atu.*W_adk));
n_e=n_gas.*P_adk;%3.32.*



end