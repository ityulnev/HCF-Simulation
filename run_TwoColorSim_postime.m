addpath(genpath('functions'))
% close all
% clear all
% initialize Parameters
beam=beam_init;
mesh=mesh_init(beam,1,1);
medium=medium_init(mesh,beam,'Neon');
pulse1=pulse_init(mesh,beam,medium,beam.Q_In,0e-15,1,0);

delays=(-60:0.5:60).*1e-15;
Ldelays=length(delays);
Ef1=zeros(Ldelays,mesh.flength);
%% Separate Fields:
pulse2=pulse_init(mesh,beam,medium,beam.Q_In2,-20e-15,2,0);
[EtF,EtS]=do_2ColorFourierSplitStep1DwithRK(mesh,beam,medium,pulse1,pulse2);

%% Separate fields Scan
for m=1:Ldelays
pulse2=pulse_init(mesh,beam,medium,beam.Q_In2,delays(m),2,0);
[EtF,EtS]=do_2ColorFourierSplitStep1DwithRK(mesh,beam,medium,pulse1,pulse2);
Ef1(m,:)=myfft(EtF,mesh);
end

%% 2DSI Trace
bb=load('2DSI_trace.mat');
Et_nir=interp1(bb.t_fs.*1e-15,bb.Et_NIR,mesh.t);
Et_nir=Et_nir./max(abs(Et_nir)).*pulse1.Epeak;
myfilterT=calc_supergaussian(mesh.t,200e-15,10,0);
Ef_nir=myfft(Et_nir.*myfilterT,mesh);
myfilterF=calc_smartfilter(mesh,Ef_nir,pulse1);
Et_nirc=myifft(2.*Ef_nir.*myfilterF,mesh);

pulse2=pulse_init(mesh,beam,medium,beam.Q_In*0.05,5e-15,1,0);%delays(m)
pulse1.Ert=Et_nirc;
[propEtF,EtFin,EtSin]=do_2Color1DwithSummedFields(mesh,beam,medium,pulse1,pulse2);
%% Initial Chirp
chirpedEt=myifft(pulse1.Erf.*exp(-1i*20.*(2*pi.*mesh.f-pulse1.w0).^2.*medium.kGVD(pulse1.pfmid)./2),mesh);
pulse1.Ert=pulse1.Epeak.*chirpedEt./max(abs(chirpedEt));

%% Summed Fields: 1D 2 Color Propagate, with Runge Kutta 

for m=1:Ldelays
pulse2=pulse_init(mesh,beam,medium,beam.Q_In,delays(m),1,0);%delays(m)

[propEtF,EtFin,EtSin]=do_2Color1DwithSummedFields(mesh,beam,medium,pulse1,pulse2);
Ef1(m,:)=myfft(propEtF,mesh);
% figure; plot(mesh.f,[abs(Ef1(m,:)).^2])
% pause(0.1) 

end
pause(0.1)

save([date,'test.mat'],'mesh','beam','medium','pulse1','delays','Ef1','Ef02','EtFin','EtSin')

