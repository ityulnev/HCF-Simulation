% close all
% clear all
% initialize Parameters
beam=beam_init;
mesh=mesh_init(beam,1);
fiber=medium_init(mesh,beam,'Neon');
pulse1=pulse_init(mesh,beam,fiber,beam.Q_In,0,1);

delays=(-60:5:60).*1e-15;
tlim=1:length(delays);
tlimMoment=1:length(delays);
%% 1D 2 Color Propagate, with Runge Kutta 
for m=1:length(delays)
pulse2=pulse_init(mesh,beam,fiber,beam.Q_In2,-400e-15,2);

[propEtF,propEtS]=do_2ColorFourierSplitStep1DwithRK(mesh,beam,fiber,pulse1,pulse2);
Ef1=myfft(propEtF,mesh);
Ef2=myfft(propEtS,mesh);
Etcomp=myifft(abs(Ef1),mesh);

tlimMoment(m)=calc_pulsedurationMoment(mesh.t,abs(Etcomp).^2);
tlim(m)=calc_fwhm(mesh.t,abs(Etcomp).^2);

calc_fwhm(mesh.t,abs(Etcomp).^2)
figure;
myplot1('intensity','no','normalize',mesh.t,Etcomp,Ecomp)
myplot1('intensity','no','no',mesh.f,Ef1,Ef2)
hold on; plot(mesh.f,[abs(IfDatFull)./max(abs(IfDatFull)).*max(abs(Efion050).^2)])

TlimStruct.(sprintf('pulse2_%d', m)) = pulse2;
TlimStruct.(sprintf('Et1_%d', m)) = propEtF;
TlimStruct.(sprintf('Et2_%d', m)) = propEtS;
TlimStruct.(sprintf('Tlim_%d', m)) = tlim(m);

calc_fwhm(mesh.t,abs(pulse1.Ert).^2)
calc_fwhm(mesh.f,abs(myfft(pulse1.Ert,mesh)).^2)


end
% save('delaycheck3_alpha_2_2mJ_1mJ.mat')

%% DPG Plots
figure;
myplot1('intensity','no','no',mesh.t.*1e15,EtnoIONSSTGVD,EtnoIONSST,EtnoION,Etall)
xlabel('time (fs)')
legend('SPM','SPM+GVD','...+SST','...+ION')
xlim([-60 60])
title('1 Color Propagation')
%%%%%%%%%%
figure;
myplot1('intensity','no','no',mesh.f.*1e-12,myfft(EtnoIONSST,mesh),myfft(EtnoION,mesh),myfft(Etall,mesh))
xlabel('frequency (THz)')
legend('SPM+GVD','...+SST','...+ION')
xlim([200 600])
title('1 Color Propagation')

Ef3=myfft(EtFsc,mesh);
myplot1('intensity','no','no',mesh.f.*1e-12,Ef3)
xlabel('frequency (THz)')
xlim([300 600])
title('Direct comparison on Log scale')

hold on; plot(mesh.f.*1e-12,[abs(IfDatFull)./max(abs(IfDatFull)).*max(abs(Ef3).^2)])
plot(mesh.f.*1e-12,[abs(IfDatFull0)./max(abs(IfDatFull0)).*max(abs(Ef3).^2)])

set(gca, 'YScale', 'log')
legend('Simulation','Measurement1','Measurement2')


%% Pulse Profile
r=(-400:1:400).*1e-6;
rmode=200e-6*0.65;
er=exp(-(((r).^2./((rmode)^2))));
er=abs(er).^2./trapz(r,abs(er).^2);

rect=zeros(1,length(r));
rect(find(r>-rmode,1):find(r>rmode,1))=1;
rect=abs(rect).^2./trapz(r,abs(rect).^2);

figure;
trapz(mesh.t,pE1)
plot(r.*1e6,[er;rect])
xlabel('r_{\perp} (µm)')
ylabel('Intensity (arb.u.)')
title('Transverse Intensity Profile')

%% DPG 2Color
Ef51=Ef1;%-52fs
Ef52=Ef2;
myplot1('intensity','no','no',mesh.f*1e-12,Ef41+Ef42)
myplot1('intensity','no','no',mesh.f*1e-12,Ef41+Ef42,Ef51+Ef52)
hold on; plot(mesh.f*1e-12,[abs(IfDatFull)./max(abs(IfDatFull)).*max(abs(Ef1+Ef2).^2)])
xlabel('frequency (THz)')
xlim([200 1000])
title('2 Color Simulation')
legend('dt=+52fs')

myplot1('intensity','no','no',mesh.f*1e-12,Ef1+Ef2)
hold on; plot(mesh.f*1e-12,[abs(IfDatFull)./max(abs(IfDatFull)).*max(abs(Ef1+Ef2).^2)])

xlabel('frequency (THz)')
legend('Simulation','Data')
xlim([300 1000])
title('Comparison with Data')


myplot1('intensity','no','no',mesh.t,pulse1.Ert,propEtF)

myplot1('intensity','no','no',mesh.t,pulse1.Ert)
set(0, 'DefaultLineLineWidth', 1);
set(0,'defaultAxesFontSize',30)
%% DPG Intensity overview

figure;
hold on
plot3(ones(1,mesh.flength),mesh.f.*1e-12,abs(Efion100).^2)
plot3(ones(1,mesh.flength).*0.9,mesh.f.*1e-12,abs(Efion090).^2)
plot3(ones(1,mesh.flength).*0.8,mesh.f.*1e-12,abs(Efion080).^2)
plot3(ones(1,mesh.flength).*0.7,mesh.f.*1e-12,abs(Efion070).^2)
plot3(ones(1,mesh.flength).*0.6,mesh.f.*1e-12,abs(Efion060).^2)
plot3(ones(1,mesh.flength).*0.5,mesh.f.*1e-12,abs(Efion050).^2)
ylim([200 600])



% trapz(mesh.f,abs(Efion080).^2)*fiber.area_hcf*fiber.Iconst

Etcomp2=myifft(abs(Efion100),mesh);
myplot1('intensity','no','no',mesh.t,Etcomp2)
calc_fwhm(mesh.t,abs(Etcomp2).^2)
S_Efion={Efion100,Efion090,Efion080,Efion070,Efion060,Efion050};



figure;
hold on
myplot1('intensity','no','no',mesh.t,pulse1.Ert)
%% DPG 2 Color Time delay Study

S_dT2Color={EtF080p40,EtS080p40,EtF080p20,EtS080p20,EtF080n00,EtS080n00,EtF080n20,EtS080n20,EtF080n40,EtS080n40};
dt=40:-20:-40;
figure;
hold on
for m=[1,3,5]
E2=pulse_init(mesh,beam,fiber,beam.Q_In2,dt(m)*1e-15,2);
plot3(ones(1,mesh.flength).*dt(m),mesh.t.*1e15,[abs(S_dT2Color{2*m-1});abs(E2.Ert);abs(S_dT2Color{2*m});S_dT2Color{2*m}])

end
ylim([-100 100])
legend('Fundamental Envelope','initial Second Harmonic','propagated Second Harminic')
%% Compression
tcomp=zeros(1,length(dt));
figure;
hold on
for m=[1,3,5]
Ef1=myfft(S_dT2Color{2*m-1},mesh);
Ef2=myfft(S_dT2Color{2*m},mesh);
Ecomp=myifft(abs(Ef1+Ef2),mesh) ;  
plot3(ones(1,mesh.flength).*dt(m),mesh.t.*1e15,abs(Ecomp).^2)    
tcomp(m)=calc_fwhm(mesh.t,abs(Ecomp).^2);    
end    
ylim([-20 20])  
figure; plot(dt,tcomp.*1e15)


figure;
hold on
for m=[1,3,5]
Ef1=myfft(S_dT2Color{2*m-1},mesh);
Ef2=myfft(S_dT2Color{2*m},mesh); 
plot(mesh.f.*1e-12,[abs(Ef2).^2])
% plot3(ones(1,mesh.flength).*dt(m),mesh.f.*1e-12,[abs(Ef1).^2;abs(Ef2).^2])    
end  
xlim([200 1200])  


% save('Test_Efion1.mat','beam','mesh','fiber','pulse1','pulse2','Efion100','Efion090','Efion100','Efion080','Efion070','Efion060','Efion050')

%% DATA FOR HCF PROPAGATION EXPERIMENT
%#############
% [x0,y0]=extract_alldata('HCFvsFund.fig');
[x1,y1,x2,y2,x3,y3,x4,y4]=extract_alldata('Fig1_2.fig');%1175:1683
% plot(x0,y0)
plot(x1,[y1;y2;y3;y4])
% [xx1,yy1]=process_realdata(x1,y1,mesh.f(1+round(mesh.flength/2):end));
[xx1,yy1]=convert_data(x3.*1e-9,y2,'back');
lbound=(find(mesh.f>xx1(end),1))
rbound=(find(mesh.f>1.5e15,1))

xb=find(xx1<1.5e15,1)

aa=interp1(xx1(xb:end),yy1(xb:end),mesh.f(lbound:rbound))
aa(isnan(aa))=0;
aanorm=aa./max(aa);

IfDatFull=[zeros(1,lbound-1),aanorm,zeros(1,mesh.flength-rbound)];
ItDat=abs(myifft(abs(IfDatFull),mesh));
figure;
hold on
plot(mesh.f.*1e-12,abs(IfDatFull))
xlim([300 600])
xlabel('frequency (THz)')
ylabel('Intensity (arb.u.)')



figure;
plot(mesh.t,abs(ItDat))
tlimDat=calc_fwhm(mesh.t,(ItDat))
% tlimDatMom=calc_pulsedurationMoment(mesh.t,ItDat);

