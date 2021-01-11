clear; clc; close all;

fnm = '2010Baja_Stress.mat';
load(fnm)
%% 
% Location 1: 20 km, SE of Parkfield, depth: 22 km
% Location 2: 30 km, NW of Parkfield, depth: 22 km
% T_Nrm: Normal Stress, tension positive
% T_Stk: Along-strike Shear Stress
% T_Dp: Along-dip Shear Stress, positive up
% t: time
% t0: date time at t = 0
% The model simulates best for the frequency range of 20 to 50 s in a sense
% of wavefrom fit for seismic observations

fs = 1/mean(diff(t));
f1 = 1/50;
f2 = 1/20;

subplot(3,1,1)
plot(t,bp_filt(T_Nrm1/1e3,f1,f2,fs),'r',t,bp_filt(T_Nrm2/1e3,f1,f2,fs),'k','linewidth',2)
title({['2010 Mw 7.1 Baja Bandpass filter at ' num2str(1/f2) ' to ' num2str(1/f1) ' s'],...
    'Normal Stress (kPa), postive tension'})
legend('20 km, SE of Parkfield','-30 km, NW of Parkfield')
set(gca,'fontsize',15)
set(gca,'Xticklabel',[])
xlim([0 400])

subplot(3,1,3)
plot(t,bp_filt(T_Stk1/1e3,f1,f2,fs),'r',t,bp_filt(T_Stk2/1e3,f1,f2,fs),'k','linewidth',2)
title('Along-strike Shear Stress (kPa), Negative shear corresponds to the righ-lateral motion')
xlabel(['Time (s), Start time: ' datestr(t0)])
legend('20 km, SE of Parkfield','-30 km, NW of Parkfield')
set(gca,'fontsize',15)
set(gca,'Xticklabel',[])
xlim([0 400])

subplot(3,1,2)
plot(t,bp_filt(T_Dp1/1e3,f1,f2,fs),'r',t,bp_filt(T_Dp2/1e3,f1,f2,fs),'k','linewidth',2)
title('Along-dip Shear Stress (kPa), positive up')
legend('20 km, SE of Parkfield','-30 km, NW of Parkfield')

set(gca,'fontsize',15)
xlim([0 400])

function wvfm_out = bp_filt(wvfm,f1,f2,fs)

bp_bp0 = designfilt('bandpassiir', 'FilterOrder',2,...
    'HalfPowerFrequency1',f1,'HalfPowerFrequency2',f2,...
    'SampleRate', fs, ...
    'DesignMethod', 'butter');
wvfm_out = filtfilt(bp_bp0,wvfm);

end