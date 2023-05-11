%% dq average model and analytical model of VSI Copyright ?2013 Boeing. All rights reserved.
clc
clear all;
close all

%% voltage loop control parameters for different cases
Kpv_VSI      = 0.0495;          % VSI voltage loop p controller           
Kiv_VSI      = 130.8997;        % VSI voltage loop i controller
Kpi_VSI      = 0.0465*270;      % VSI current loop p controller
Kii_VSI      = 121.0095*270;    % VSI current loop i controller

%% opertions for Bode plot display
Bode_O=bodeoptions;
Bode_O.XLabel.FontSize=14;
Bode_O.YLabel.FontSize=14;
Bode_O.TickLabel.FontSize=14;
Bode_O.Title.FontSize=14;
Bode_O.title.String=' ';
Bode_O.Grid='on';
Bode_O.XLim={[40 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='on';
%% Matlab settings, in case encoding doesn't work
% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
%% converter power stage parameters
L=1.0e-3;
RL=0.13;
R=12;
C=31.5e-6;
Cdc=150e-6;
f=400;
w=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;

Vdc=270; Vdcref=270;
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsdq=[sqrt(3/2)*Vsm; 0];
Vdref=Vsdq(1); Vqref=Vsdq(2);

%% model the PWM delay both for both average model and analytical model
s=tf([1,0],[0,1]);
I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
Tdelay = 1.5/fsw;
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
sys = linearize('Zdc_VSI_v_loop',0.5);
%% Zin model linearization
H=ss(sys);
Vavg=H(1,1);
Ilavg=H(2,1);
ZVSI_dc=Vavg/Ilavg;

%% impedance comparison
fig=figure(1)
set(fig, 'Position', [50, 10, 1000, 800]);
bode(ZVSI_dc,Bode_O);
legend('ZVSIdc\_avg\_sim')
Bode_Darklines(3)