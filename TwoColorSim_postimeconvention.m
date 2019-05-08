addpath(genpath('functions'))
close all
clear all
% initialize Parameters
beam=beam_init;
mesh=mesh_init(beam,1);
fiber=medium_init(mesh,beam,'Neon');
pulse1=pulse_init(mesh,beam,fiber,beam.Q_In,0,1);

delays=(-20:0.1:20).*1e-15;
Ef1=zeros(length(delays),mesh.flength);
Ef2=zeros(length(delays),mesh.flength);

%% 1D 2 Color Propagate, with Runge Kutta 
for m=1:1:length(delays)
pulse2=pulse_init(mesh,beam,fiber,beam.Q_In2,delays(m),2);

[propEtF,propEtS]=do_2ColorFourierSplitStep1DwithRK(mesh,beam,fiber,pulse1,pulse2);
% Ef1=myfft(propEtF,mesh);
% Ef2=myfft(propEtS,mesh);

Ef1(m,:)=myfft(propEtF,mesh);
Ef2(m,:)=myfft(propEtS,mesh);

% plot(mesh.f,[Ef1(m,:);abs(Ef1(m,:));Ef2(m,:);abs(Ef2(m,:))])
% pause(0.1)
% Etcomp=myifft(abs(Ef1),mesh);

end
save([date,'_GVDSPMdelaycheck.mat'],'mesh','beam','fiber','pulse1','delays','Ef1','Ef2')
