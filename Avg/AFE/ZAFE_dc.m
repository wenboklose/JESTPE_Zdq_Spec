%% dq average model and analytical model of AFE Copyright ?2013 Boeing. All rights reserved.
clc
clear all;
close all
%% voltage and current control loop as well as PLL parameters for AFE which is fixed for all the cases
Kpv_AFE      = 0.0462;          % AFE voltage loop p controller
Kiv_AFE      = 4.5815;          % AFE voltage loop i controller
Kpi_AFE      = 0.0116*270;      % AFE current loop p controller
Kii_AFE      = 46.5421*270;     % AFE current loop i controller
FWPI_a1      = 8.9207;          % AFE Phase-Locked Loop pi controller
FWPI_a2      = -8.7225;         % AFE Phase-Locked Loop pi controller
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
Vdcref=270;
Cdc=105e-6; 
R=96;
RCdc=0.049;
P=Vdcref^2/R;
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd_s=Vsdq(1);
Iq_ref=0;
L=471e-6;
RL=110e-3;
f=400;
omega=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
% model the PWM delay both for both average model and analytical model
s=tf([1,0],[0,1]);
I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
Tdelay = 1.5/fsw;
Gsvm = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);
Gsvm = [Gsvm 0*Gsvm; -0*Gsvm Gsvm];
Gsvm = tf(JF_DQFromABC(Gsvm(1,1),omega));
Gsvm = [Gsvm(1,1) 0*Gsvm(1,2); -0*Gsvm(1,2) Gsvm(2,2)];
Gdel = Gsvm;
%% model the signal conditioning filter for both average model and analytical model
C1=75e-12;C2=39e-12;R1=18e3;R2=18e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));
tf_filter_dq = [tf_filter_dq(1,1) 0*tf_filter_dq(1,2); -0*tf_filter_dq(1,2) tf_filter_dq(2,2)];
Ki = tf_filter_dq;
Kv = tf_filter_dq;

%% assign control and PLL parameters for models

kpi=Kpi_AFE/270;
kii=Kii_AFE/270;
kpv=Kpv_AFE;
kiv=Kiv_AFE;
tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],Tsw);
tf_pll=d2c(tf_pll_z);

%% dq average model linearization:
sys = linearize('ZAFE_dc_v_loop',0.5);
%% Zin model linearization
H=ss(sys);
Vsvg1=H(2)
Isavg1=H(1)
Zdc_AFE=Vsvg1/Isavg1;

fig=figure(2)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(Zdc_AFE,Bode_O)
legend('ZAFEdc\_avg\_sim')
Bode_Darklines(3)