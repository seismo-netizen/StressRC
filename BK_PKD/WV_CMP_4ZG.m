clear; clc; close all;

fs = 1/0.02;
f1 = 1/50;
f2 = 1/20;

tarnm = 'RC64.BK.PKD'; step_v = 0.05; step_a = 0.02;
%tarnm = 'RC71.BK.PKD'; step_v = 0.25; step_a = 0.1;

%% Velocity
figure
subplot(1,2,1)
sacfnm = [tarnm '.HHZ.00.sac'];
[hdr,data] = load_sac(sacfnm); w = tukeywin(length(data),0.05);
data = bp_filt(data.*w,f1,f2,1/hdr.delta);
t2 = 1:hdr.npts; t2 = t2 * hdr.delta + hdr.b;
plot(t2,data,'-k','linewidth',2); hold on;

sacfnm = [tarnm '.HHN.00.sac'];
[hdr,data] = load_sac(sacfnm);
data = bp_filt(data.*w,f1,f2,1/hdr.delta);
t2 = 1:hdr.npts; t2 = t2 * hdr.delta + hdr.b;
plot(t2,data-step_v,'-k','linewidth',2)

sacfnm = [tarnm '.HHE.00.sac'];
[hdr,data] = load_sac(sacfnm);
data = bp_filt(data.*w,f1,f2,1/hdr.delta);
t2 = 1:hdr.npts; t2 = t2 * hdr.delta + hdr.b;
plot(t2,data-step_v*2,'-k','linewidth',2)


load([tarnm '.mat']);
dt = mean(diff(t)); w = tukeywin(length(t)-1,0.05);
Vz = diff(rcv_z)/dt; Vz = bp_filt(Vz.*w',f1,f2,1/dt);
Ve = diff(rcv_x)/dt; Ve = bp_filt(Ve.*w',f1,f2,1/dt);
Vn = diff(rcv_y)/dt; Vn = bp_filt(Vn.*w',f1,f2,1/dt);

plot(t(1:end-1), Vz*100,'-r','linewidth',2);
plot(t(1:end-1), Vn*100-step_v,'-r','linewidth',2);
plot(t(1:end-1), Ve*100-step_v*2,'-r','linewidth',2); 
xlim([0 300]); ylim([-step_v*2.5 step_v*0.5])

text(-40,0,'HHZ','Fontsize',15)
text(-40,-step_v,'HHN','Fontsize',15)
text(-40,-step_v*2,'HHE','Fontsize',15)

plot([200 200],[-0.25*step_v -0.75*step_v],'k','linewidth',2)
text(210,-0.5*step_v,[num2str(step_v/2) 'cm/s'],'Fontsize',15)


set(gca,'Yticklabel',[])
xlabel('Time (s)')
title([tarnm ' Velocity (cm/s)'])
set(gca,'fontsize',15)

%% Acceleration
subplot(1,2,2)
sacfnm = [tarnm '.HNZ.00.sac'];
[hdr,data] = load_sac(sacfnm); w = tukeywin(length(data),0.05);
data = bp_filt(data.*w,f1,f2,1/hdr.delta);
t2 = 1:hdr.npts; t2 = t2 * hdr.delta + hdr.b;
plot(t2,data,'-k','linewidth',2); hold on;

sacfnm = [tarnm '.HNN.00.sac'];
[hdr,data] = load_sac(sacfnm);
data = bp_filt(data.*w,f1,f2,1/hdr.delta);
t2 = 1:hdr.npts; t2 = t2 * hdr.delta + hdr.b;
plot(t2,data-step_a,'-k','linewidth',2)

sacfnm = [tarnm '.HNE.00.sac'];
[hdr,data] = load_sac(sacfnm);
data = bp_filt(data.*w,f1,f2,1/hdr.delta);
t2 = 1:hdr.npts; t2 = t2 * hdr.delta + hdr.b;
plot(t2,data-step_a*2,'-k','linewidth',2)

Az = diff(Vz)/dt; 
Ae = diff(Ve)/dt; 
An = diff(Vn)/dt;

plot(t(1:end-2), Az*100,'-r','linewidth',2); hold on;
plot(t(1:end-2), An*100-step_a,'-r','linewidth',2); hold on;
plot(t(1:end-2), Ae*100-step_a*2,'-r','linewidth',2); hold on;

text(-40,0,'HNZ','Fontsize',15)
text(-40,-step_a,'HNN','Fontsize',15)
text(-40,-step_a*2,'HNE','Fontsize',15)

plot([200 200],[-0.25*step_a -0.75*step_a],'k','linewidth',2)
text(210,-0.5*step_a,[num2str(step_a/2) 'cm/s^2'],'Fontsize',15)

xlim([0 300]); ylim([-step_a*2.5 step_a*0.5])
set(gca,'Yticklabel',[])
xlabel('Time (s)')
title([tarnm ' Acceleration (cm/s^2)'])
set(gca,'fontsize',15)


function wvfm_out = bp_filt(wvfm,f1,f2,fs)

bp_bp0 = designfilt('bandpassiir', 'FilterOrder',4,...
    'HalfPowerFrequency1',f1,'HalfPowerFrequency2',f2,...
    'SampleRate', fs, ...
    'DesignMethod', 'butter');
wvfm_out = filtfilt(bp_bp0,wvfm);

end