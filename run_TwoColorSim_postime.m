addpath(genpath('functions'))
% close all
% clear all
% initialize Parameters
beam=beam_init;
mesh=mesh_init(beam,1,1);
% medium=medium_init(mesh,beam,'Neon');
medium=gen_medium_init(mesh,beam,'Neon');
% pulse1=pulse_init(mesh,beam,medium,beam.Q_In,0e-15,1,0);

% delays=(-40:5:40).*1e-15;
% Ldelays=length(delays);
% Ef1=zeros(Ldelays,mesh.flength);
% Ef2=zeros(Ldelays,mesh.flength);

%% Separate Fields:
r_mode=0.65*200e-6;
beam_area=pi*r_mode^2;
t_pulse=35e-15;
t_fwhm=t_pulse.*(sqrt(log(2)/2));
Ipeak=1e-3*0.94/(t_fwhm*beam_area/2);
pulse1=general_pulse_init(mesh,800e-9,t_pulse,r_mode,Ipeak,medium.Iconst,0,0);
% for m=1:Ldelays
pulse2=general_pulse_init(mesh,800e-9/3,t_pulse,r_mode,Ipeak*0.1,medium.Iconst,0e-15,0);
pulse2m=general_pulse_init(mesh,800e-9/3,t_pulse,r_mode,Ipeak*0.1,medium.Iconst,20e-15,0);
pulse2p=general_pulse_init(mesh,800e-9/3,t_pulse,r_mode,Ipeak*0.1,medium.Iconst,-15e-15,0);
% figure; plot(mesh.t,abs([pulse1.Ert;pulse2.Ert]).^2)
[Ef11,Ef21]=do_2ColorFourierSplitStep1DwithRK(mesh,beam,medium,pulse1,pulse2m);
% save('xpm_phasestudy.mat','myS')
Et21=myifft(cmpns_tshift(Ef21,mesh),mesh);

% myf0=1200e12;
myf0=calc_centerofmass(mesh.f,abs(Ef21).^2','cartesian')
index_myf0=find(mesh.f>myf0,1);
mat1=material_init(mesh,'FusedSilica');
k2=mat1.kGVD(index_myf0);
gvd=-k2*2*pi.*(mesh.f-mesh.f(index_myf0)).^2.*2.5e-4./2;
phase_gvd=exp(1i.*gvd);
phase_n=exp(-1i.*(mat1.k-mat1.kGV(index_myf0)*2*pi.*(mesh.f-myf0)).*0.35e-4); 
phase_n=handle_NaNInf(phase_n);
phase_n(abs(phase_n)>1000)=0;

mat2=material_init(mesh,'Air');

myfilt=calc_supergaussian(mesh.f,300e12,10,myf0);%pulse2m.f0
figure; plot(mesh.f.*1e-12,[abs(Ef21+Ef11).^2./max(abs(Ef21+Ef11).^2);myfilt;mat1.n-1;(mat2.n-1).*1000]); xlim([100,2000])
Et21_compr_orig=myifft(cmpns_tshift(Ef21,mesh).*phase_n,mesh);
Et21_compr=myifft((cmpns_tshift(Ef21+Ef11,mesh)).*phase_n.*myfilt,mesh);
Et21_compr_air=myifft((cmpns_tshift(Ef21+Ef11,mesh)).*myfilt.*exp(-1i.*mat2.k.*0.06),mesh);
hold on; plot(mesh.f.*1e-12,[abs(myfft(Et21_compr,mesh)).^2./max(abs(Ef21+Ef11).^2)])
calc_fwhm(mesh.t,abs(Et21_compr_orig).^2)
calc_fwhm(mesh.t,abs(Et21_compr).^2)
calc_fwhm(mesh.t,abs(Et21_compr_air).^2)
figure; plot(mesh.t.*1e15,[abs(Et21_compr_orig).^2;abs(Et21_compr).^2;abs(Et21_compr_air).^2])

%% convert to wavelength


%% plot
figure; 
% subplot(2,1,1)
plot(mesh.t.*1e15,pulse1.Iconst.*1e-4.*[abs(Et21).^2;abs(Et21_compr).^2],'LineWidth',1)
xlim([-50,50]); xlabel('time (fs)');
ylabel('I (W/cm$^2$)'); ylim([0,12e13])
legend('uncompressed','compressed')
% my_figure_settings('xpm_It_-15fs_FS25e-5',1)
subplot(2,1,2)
yyaxis left
plot(mesh.f.*1e-12,[abs(Ef27).^2])
yyaxis right
phase1=unwrap(angle(cmpns_tshift(Ef27,mesh)));
plot(mesh.f.*1e-12,[(phase1-phase1(index_myf0));-gvd])
xlim([600,1600])

figure; plot(mesh.t.*1e15,[real(Et21_compr);abs(Et21_compr)])




%% Separate fields Scan
for m=1:Ldelays
pulse2=pulse_init(mesh,beam,medium,beam.Q_In2,delays(m),2,0);
[EtF,EtS]=do_2ColorFourierSplitStep1DwithRK(mesh,beam,medium,pulse1,pulse2);
Ef1(m,:)=myfft(EtF,mesh)+myfft(EtS,mesh);
% Ef2(m,:)=myfft(EtS,mesh);
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
pulse2=pulse_init(mesh,beam,medium,beam.Q_In,-20e-15,1,0);%delays(m)
[propEtF,EtFin,EtSin]=do_2Color1DwithSummedFields(mesh,beam,medium,pulse1,pulse2);
Ef1(m,:)=myfft(propEtF,mesh);
% figure; plot(mesh.f,[abs(Ef1(m,:)).^2])
% pause(0.1) 
end
pause(0.1)

% save([date,'test.mat'],'mesh','beam','medium','pulse1','delays','Ef1','Ef02','EtFin','EtSin')

