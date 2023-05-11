%% dq average model and analytical model of AFE Copyright ?2013 Boeing. All rights reserved.
clc
clear all;
close all
%% voltage and current control loop as well as PLL parameters for AFE which is fixed for all the cases

%% opertions for Bode plot display
Bode_O=bodeoptions;
Bode_O.XLabel.FontSize=14;
Bode_O.YLabel.FontSize=14;
Bode_O.TickLabel.FontSize=14;
Bode_O.Title.FontSize=14;
Bode_O.title.String=' ';
Bode_O.Grid='on';
Bode_O.XLim={[1 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='on';
%% Matlab settings, in case encoding doesn't work
% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
%% filter parameters:
Lf=120e-6;
RLf=150e-3;
Cf=10e-6;
Rcf=50e-3;
%% converter power stage parameters
Vdcref=270;%320;%
Cdc=105e-6; 
R=96;%40;%50;%60;%70;%80;%
RCdc=0.049;
P=Vdcref^2/R;
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd_s=Vsdq(1);
Iq_ref=0;%-7;%
L=470e-6;
RL=100e-3;
f=60;%400;
omega=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
% model the PWM delay both for both average model and analytical model
s=tf([1,0],[0,1]);
I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
Tdelay = 1.5/fsw;%1.5/fsw;
Gsvm = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);
Gsvm = [Gsvm 0*Gsvm; -0*Gsvm Gsvm];
Gsvm = tf(JF_DQFromABC(Gsvm(1,1),omega));
Gsvm = [Gsvm(1,1) 0*Gsvm(1,2); -0*Gsvm(1,2) Gsvm(2,2)];
Gdel = Gsvm;
%% model the signal conditioning filter for both average model and analytical model
C1=75e-12;C2=39e-12;R1=18e3;R2=18e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));
% tf_filter_dq = [tf_filter_dq(1,1) 0*tf_filter_dq(1,2); -0*tf_filter_dq(1,2) tf_filter_dq(2,2)];
Ki = tf_filter_dq;
Kv = tf_filter_dq;

%% assign control and PLL parameters for models
fv=100; 
fi=1000;
Lboost_con=500e-6;
RLboost_con=100e-3;
Cdc_con=100e-6;
Rdc_con=96;
kpi=2*pi*fi*Lboost_con/Vdcref;
kii=10*2*pi*fi*RLboost_con/Vdcref;
kpv=2*pi*fv*Cdc_con;
kiv=0.5*2*pi*fv/Rdc_con;

DEF_pll=10;%40;%                     		%control loop for PLL (Hz)
DEF_pll_damp=0.707;% 0.4;%                   %damping factor for the PLL controller
DEF_Vin=57.5;                           %input phase to neutral rms voltage
DEF_Tsw=0.00005;						%switching period
FWPI_a1=2*DEF_pll_damp*DEF_pll*6.28318530717959*0.57735026918963/DEF_Vin;					%digital PI parameters for frequency loop
FWPI_a2=-(2*DEF_pll_damp*DEF_pll*6.28318530717959-DEF_pll*DEF_pll*39.47841760435743*DEF_Tsw)*0.57735026918963/DEF_Vin;
tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],DEF_Tsw);
tf_pll=d2c(tf_pll_z);

%% dq average model linearization:
% sys = linearize('ZAFE_ac_v_loop',0.5);
% %% Zin model linearization
% H=ss(sys);
% Vavg1=[H(1,1) H(1,2);
%     H(2,1) H(2,2);];
% Ilavg1=[H(3,1) H(3,2);
%     H(4,1) H(4,2);];
% Zac_AFE=Vavg1/Ilavg1;

%% dq average i decoupled model linearization:
sys = linearize('ZAFE_ac_v_loop_i_dec',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];
Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];
Zac_AFE_i_dec=Vavg1/Ilavg1;
% %% dq average i decoupled model linearization:
% sys = linearize('ZAFE_ac_v_loop_i_dec_LC_filter',0.5);
% %% Zin model linearization
% H=ss(sys);
% Vavg1=[H(1,1) H(1,2);
%     H(2,1) H(2,2);];
% Ilavg1=[H(3,1) H(3,2);
%     H(4,1) H(4,2);];
% Zac_AFE_i_dec_filter=Vavg1/Ilavg1;
% 
Yac_AFE_i_dec=inv(Zac_AFE_i_dec);
fig=figure(1)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(Yac_AFE_i_dec(1,1),Yac_AFE_i_dec(1,2),Yac_AFE_i_dec(2,1),Yac_AFE_i_dec(2,2),Bode_O)
legend('ZAFEdq\_avg\_sim')
Bode_Darklines(3)




%% model for EMI filter:
Ldm = 120e-6;
RLdm = 1*150e-3;
Cdm = 10e-6;
RCdm = 2*0.004/(omega*Cdm);
s = tf([1 0],[0 1]);
ZLdm = [Ldm*s + RLdm -Ldm*60*2*pi;Ldm*60*2*pi Ldm*s + RLdm];

YCdm = [Cdm*s -Cdm*60*2*pi;Cdm*60*2*pi Cdm*s];
ZCdm = [RCdm 0;0 RCdm]+1/YCdm;
Yac_AFE_i_dec_filter=inv(ZLdm+inv(1/ZCdm+Yac_AFE_i_dec));
Yfilter=inv(ZLdm+ZCdm);
Zfilter=ZLdm+ZCdm;
Zfilter_O=1/(1/ZLdm+1/ZCdm);
% Yac_AFE_i_dec_filter=inv(Zac_AFE_i_dec_filter);

fig=figure(2)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(Yac_AFE_i_dec_filter(1,1),Yac_AFE_i_dec_filter(1,2),Yac_AFE_i_dec_filter(2,1),Yac_AFE_i_dec_filter(2,2),Bode_O)
legend('ZAFEdq\_avg\_sim')
Bode_Darklines(3)


% fig=figure(21)
% set(fig, 'Position', [50, 10, 1000, 800]);
% bode(Zfilter_O(1,1),Zfilter_O(1,2),Bode_O)
% legend('ZAFEdq\_avg\_sim')
% Bode_Darklines(3)
% 
% % load ZAFE_270_60_pll_10_20_10k.mat
% % Zac_AFE_meas=ZAFE_320_60_20_10k;
% load ZAFE_270_60_pll_40_20_10k.mat
% Zac_AFE_meas=ZAFE_270_60_pll_40_20_10k;
% Yac_AFE_meas=inv(Zac_AFE_meas);
% fig=figure(3)
% set(fig, 'Position', [50, 10, 1000, 800]);
% bode(1/Zac_AFE_i_dec,1/Zac_AFE_meas,Bode_O)
% Bode_Darklines(3)
% 
% 
% fig=figure(41)
% set(fig, 'Position', [50, 10, 1000, 800]);
% bode(Zac_AFE_i_dec(1,1),Zac_AFE_meas(1,1),Bode_O)
% Bode_Darklines(3)
% 
% fig=figure(42)
% set(fig, 'Position', [50, 10, 1000, 800]);
% bode(Zac_AFE_i_dec(1,2),Zac_AFE_meas(1,2),Bode_O)
% Bode_Darklines(3)
% 
% fig=figure(43)
% set(fig, 'Position', [50, 10, 1000, 800]);
% bode(Zac_AFE_i_dec(2,1),Zac_AFE_meas(2,1),Bode_O)
% Bode_Darklines(3)
% 
% fig=figure(44)
% set(fig, 'Position', [50, 10, 1000, 800]);
% bode(Zac_AFE_i_dec(2,2),Zac_AFE_meas(2,2),Bode_O)
% Bode_Darklines(3)
% 
% 
% 
% fig=figure(51)
% set(fig, 'Position', [50, 10, 1000, 800]);
% bode(Yac_AFE_i_dec(1,1),Yac_AFE_meas(1,1),Bode_O)
% Bode_Darklines(3)
% 
% fig=figure(52)
% set(fig, 'Position', [50, 10, 1000, 800]);
% bode(Yac_AFE_i_dec(1,2),Yac_AFE_meas(1,2),Bode_O)
% Bode_Darklines(3)
% 
% fig=figure(53)
% set(fig, 'Position', [50, 10, 1000, 800]);
% bode(Yac_AFE_i_dec(2,1),Yac_AFE_meas(2,1),Bode_O)
% Bode_Darklines(3)
% 
% fig=figure(54)
% set(fig, 'Position', [50, 10, 1000, 800]);
% bode(Yac_AFE_i_dec(2,2),Yac_AFE_meas(2,2),Bode_O)
% Bode_Darklines(3)