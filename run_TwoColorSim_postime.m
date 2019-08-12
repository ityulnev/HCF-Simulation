addpath(genpath('functions'))
% close all
% clear all
% initialize Parameters
beam=beam_init;
mesh=mesh_init(beam,1,1);
medium=medium_init(mesh,beam,'Neon');
pulse1=pulse_init(mesh,beam,medium,beam.Q_In,0e-15,1,0);
% pulse3=pulse_init(mesh,beam,fiber,beam.Q_In,10e-15,1,0);
%% Separate Fields:
pulse2=pulse_init(mesh,beam,medium,beam.Q_In*0.1,200e-15,2,0);
[EtF,EtS]=do_2ColorFourierSplitStep1DwithRK(mesh,beam,medium,pulse1,pulse2);

%% Summed Fields: 1D 2 Color Propagate, with Runge Kutta 
delays=(-30 :0.2:30).*1e-15;
Ef1=zeros(length(delays),mesh.flength);

for m=1:1:length(delays)
pulse2=pulse_init(mesh,beam,medium,beam.Q_In,20e-15,1,0);%delays(m)

[propEtF]=do_2Color1DwithSummedFields(mesh,beam,medium,pulse1,pulse2);
Ef02(m,:)=pulse2.Erf;%history of second input pulse
Ef1(m,:)=myfft(propEtF,mesh);
% figure; plot(mesh.f,[abs(Ef1(m,:)).^2])
% pause(0.1)

end

save([date,'_GVDSPMdelaycheck.mat'],'mesh','beam','fiber','pulse1','delays','Ef1','Ef2','Ef02')

