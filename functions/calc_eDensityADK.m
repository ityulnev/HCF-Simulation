%Calculates W_adk for a given E field and a given Gas
function [n_e,Eg]=calc_eDensityADK(E_laser,mesh,fiber)
E_laser=abs(E_laser);
%Quantum Numbers
l=1;
m=0;
f=(2*l+1)*factorial(l+abs(m))/(2^abs(m)*factorial(abs(m))*factorial(l-m));
Z=1;                                                                       %First ionization
%Gastype
switch fiber.gas
    case 'Neon'
            Eg=21.565*const.e;                                             %[J] from http://www.periodensystem.info/elemente/neon/
            n_gas=2.686e25*fiber.pressure;                                                %1/m^3 
    case 'Argon'
            Eg=15.76*const.e;                                              %[J] from http://www.periodensystem.info
            n_gas=2.7e25*fiber.pressure;                                                  %1/m^3   
    case 'Xenon'
            Eg=12.13*const.e;                                              %[J] from http://www.periodensystem.info
            n_gas=2.4e25*fiber.pressure;                                                  %1/m^3  
end 
%% Conversion into Atomic Units
rbohr=4*pi*const.eps0*const.hbar^2/(const.m_e*const.e^2);   
E_hartree=const.m_e*const.e^4/(4*const.eps0^2*const.h^2);                  %Hartree Energy [J]
Efield_hartree=E_hartree/(rbohr*const.e);                                  %Efield from Hartree Energy [V/m]
Eg_au=Eg/E_hartree;                                                        %Energy into atomic units a.u.
E_laser_au=E_laser./Efield_hartree;                                        %Field into atomic untis a.u.  
dt_atu=mesh.dt./(const.hbar/E_hartree);                             %Time into atomic time units a.t.u.

%% Calc ionization density n_e with ADK model
n_star=Z*(1/sqrt(2*Eg_au));                                                %Effective principal quantum number
E_0_au=(2*Eg_au)^(3/2);
fraction=E_0_au./E_laser_au;
exponent=exp(-(2/3).*fraction);
w_factor=(E_laser_au).^(1.5-2*n_star);
w_factor(isinf(w_factor))=0;
prod=w_factor.*exponent;

W_adk=f/(8*pi*n_star).*((4*exp(1)*E_0_au./(n_star)).^(2*n_star)).*sqrt(3./(pi*E_0_au*(2*Eg_au))).*prod;  %Ionization rate
P_adk=1-exp(-dt_atu.*W_adk);                                               %Ioniz. propability
n_e=cumsum(n_gas.*(P_adk));                                                %Electron density

end