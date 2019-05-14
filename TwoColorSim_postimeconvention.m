addpath(genpath('functions'))
% close all
% clear all
% initialize Parameters
beam=beam_init;
mesh=mesh_init(beam,1);
fiber=medium_init(mesh,beam,'Neon');
pulse1=pulse_init(mesh,beam,fiber,beam.Q_In,0,1);

delays=(-15:0.1:15).*1e-15;
Ef1=zeros(length(delays),mesh.flength);

%% 1D 2 Color Propagate, with Runge Kutta 
for m=1:1:length(delays)
pulse2=pulse_init(mesh,beam,fiber,beam.Q_In2,delays(m),1);%delays(m)

[propEtF]=do_2Color1DwithSummedFields(mesh,beam,fiber,pulse1,pulse2);
Ef02(m,:)=pulse2.Erf;%history of second input pulse
Ef1(m,:)=myfft(propEtF,mesh);
% figure; plot(mesh.f,[abs(Ef1(m,:)).^2])
% pause(0.1)

end
save([date,'_GVDSPMdelaycheck.mat'],'mesh','beam','fiber','pulse1','delays','Ef1','Ef2','Ef02')

