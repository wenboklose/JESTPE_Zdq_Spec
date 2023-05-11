%% dq average model and analytical model of VSI Copyright ?2013 Boeing. All rights reserved.
clc
clear all;
close all

%% voltage loop control parameters for different cases

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
%% converter power stage parameters
L=970e-6;
RL=0.11;
R=12;
RLoad = 215;
C=31.8e-6;
CR=1;
Cdc=150e-6;
f=60;%400;%
w=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;

Vdc=270; Vdcref=270;
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsdq=[sqrt(3/2)*Vsm; 0];
Vdref=Vsdq(1); Vqref=Vsdq(2);

%% calculate control parameters
fv=150;%150;%1500;
fi=3500;%1000;
Lcon=1e-3;
RLcon = 150e-3;
Ccon=30.8e-6;
Rcon=12;
Kpi_VSI=.1*2*pi*fi*Lcon*0.5;
Kii_VSI=0*2*pi*fi*(RLcon)*10;
Kpv_VSI=2*pi*fv*Ccon;
Kiv_VSI=2*pi*fv/Rcon;
K_dec =1*0.95;

%% model the PWM delay both for both average model and analytical model
s=tf([1,0],[0,1]);
I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
Tdelay = 1.5/fsw;%1.5/fsw;
Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);
Gsvm = tf(JF_DQFromABC(Gsvmdd,w));
Gsvm = [Gsvm(1,1) 0*Gsvm(1,2); -0*Gsvm(1,2) Gsvm(2,2)];
Gdel = Gsvm;
%% model the signal conditioning filter for both average model and analytical model
C1=75e-12;C2=39e-12;R1=18e3;R2=18e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,w));
tf_filter_dq = [tf_filter_dq(1,1) 0*tf_filter_dq(1,2); -0*tf_filter_dq(1,2) tf_filter_dq(2,2)];
Ki = tf_filter_dq;
Kv = tf_filter_dq;

%% assign control parameters for models

kpi=Kpi_VSI/270;
kii=Kii_VSI/270;
kpv=Kpv_VSI;
kiv=Kiv_VSI;

%% dq average model linearization:
sys = linearize('Zac_VSI_v_loop',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];

Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];

Isavg1=[H(5,1) H(5,2);
    H(6,1) H(6,2);];
ZVSI_ac=Vavg1/Isavg1;
Dd=0.3674;
Dq=0.0012;
%% dq average v decouple model linearization:
sys = linearize('VSI_o_loop_Z_ac',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];

Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];

Isavg1=[H(5,1) H(5,2);
    H(6,1) H(6,2);];
ZVSI_ac_o=Vavg1/Isavg1;
%% dq average v decouple model linearization:
sys = linearize('VSI_o_loop_Z_ac_Rdamp_C',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];

Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];

Isavg1=[H(5,1) H(5,2);
    H(6,1) H(6,2);];
ZVSI_ac_o_Rdamp =Vavg1/Isavg1;
% dq average i v decouple model linearization:
sys = linearize('Zac_VSI_v_loop_i_v_dec',0.5);
%% Zin model linearization
H=ss(sys);
Vavg1=[H(1,1) H(1,2);
    H(2,1) H(2,2);];

Ilavg1=[H(3,1) H(3,2);
    H(4,1) H(4,2);];

Isavg1=[H(5,1) H(5,2);
    H(6,1) H(6,2);];
ZVSI_ac_i_v_dec=Vavg1/Isavg1;

% Dd=0.3674;
% Dq=0.0012;
% sys = linearize('VSI_o_loop_Z_ac',0.5);
% 
% %% Zin model linearization
% H=ss(sys);
% Vavg1=[H(1,1) H(1,2);
%     H(2,1) H(2,2);];
% 
% Ilavg1=[H(3,1) H(3,2);
%     H(4,1) H(4,2);];
% 
% Isavg1=[H(5,1) H(5,2);
%     H(6,1) H(6,2);];
% ZVSI_ac_o=Vavg1/Isavg1;

% impedance comparison
%% open loop without damping
fig=figure(1)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(ZVSI_ac_o(1,1),ZVSI_ac_o(1,2),Bode_O);
Bode_Darklines(3)
%% open loop with damping
fig=figure(2)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(ZVSI_ac_o_Rdamp(1,1),ZVSI_ac_o_Rdamp(1,2),Bode_O);
Bode_Darklines(3)
%% closed-loop with damping and decoupling terms
fig=figure(3)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(ZVSI_ac_i_v_dec(1,1),ZVSI_ac_i_v_dec(1,2),Bode_O);
Bode_Darklines(3)

AZ = 1/sqrt(2)*[1 1i;1 -1i];
ZVSI_pn = AZ*ZVSI_ac_i_v_dec*inv(AZ);

fig=figure(31)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(ZVSI_pn(1,1),ZVSI_pn(2,2),Bode_O);
Bode_Darklines(3)

%% ZVSI with cable (inductive impedance, non decoupled)
%% impedance of the cable (1km):
% low voltage: R: 0.642Ohm/km X: 0.083Ohm/km
omega=2*pi*f;
meter=1;  % of 1km
RL = meter*0.642;
L = meter*0.083/(omega);
% medium voltage: R: 0.161Ohm/km X: 0.190Ohm/km
RL = meter*0.161;
L = meter*0.19/omega;
% high voltage: R: 0.06Ohm/km X: 0.191Ohm/km
RL = meter*0.06;
L = meter*0.191/omega
% C = 10e-6;
% RC = 2*0.000004/(omega*C);
s = tf([1 0],[0 1]);
ZLdq = [L*s + RL -L*omega;L*omega L*s + RL];

ZVSI_L_dq = ZVSI_ac_i_v_dec + ss(ZLdq);
fig=figure(4)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(ZVSI_L_dq(1,1),ZVSI_L_dq(1,2),Bode_O);
Bode_Darklines(3)

n=1e2;
f=logspace(-1,4,n);
w=f*2*pi;

Zdd = freqresp(ZVSI_L_dq(1,1),w);
Zdq = freqresp(ZVSI_L_dq(1,2),w);
Zqd = freqresp(ZVSI_L_dq(2,1),w);
Zqq = freqresp(ZVSI_L_dq(2,2),w);

Zdd_frd = frd(Zdd,w);
Zdq_frd = frd(Zdq,w);
Zqd_frd = frd(Zqd,w);
Zqq_frd = frd(Zqq,w);
ZVSI_L_dq_frd = [Zdd_frd Zdq_frd;Zqd_frd Zqq_frd];

AZ = 1/sqrt(2)*[1 1i;1 -1i];
ZVSI_pn = AZ*ZVSI_L_dq_frd*inv(AZ);

fig=figure(41)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(ZVSI_pn(1,1),ZVSI_pn(2,1),Bode_O);
Bode_Darklines(3)

    
load ZVSI_60_20_10k_exp.mat
ZVSI_60_20_10k_meas=ZVSI_60_20_10k;
fig=figure(5)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(ZVSI_ac_i_v_dec,ZVSI_60_20_10k,Bode_O);
Bode_Darklines(3)

% fig=figure(3)
% set(fig, 'Position', [50, 10, 1000, 800]);
% nyquist(ZVSI_ac_i_v_dec,[0.1*2*pi:10*2*pi:10000*2*pi]);
% Bode_Darklines(3)
% 
% RGA_ZVSI = ZVSI_ac_i_v_dec*inv(ZVSI_ac_i_v_dec)';
% fig=figure(4)
% bode(RGA_ZVSI(1,1),RGA_ZVSI(1,2),Bode_O);
% Bode_Darklines(3)