clc
clear all;
close all
%% voltage and current control loop as well as PLL parameters for AFE which is fixed for all the cases
Kpv_AFE      = 2*0.0462;          % AFE voltage loop p controller
Kiv_AFE      = 2*4.5815;          % AFE voltage loop i controller
Kpi_AFE      = 1*0.0116*270;      % AFE current loop p controller
Kii_AFE      = 1*46.5421*270;     % AFE current loop i controller
FWPI_a1      = 8.9207/10;          % AFE Phase-Locked Loop pi controller
FWPI_a2      = -8.7225/10;         % AFE Phase-Locked Loop pi controller
%% Bode plot options
Bode_O=bodeoptions;
Bode_O.XLabel.FontSize=18;
Bode_O.YLabel.FontSize=18;
Bode_O.TickLabel.FontSize=18;
Bode_O.Title.FontSize=18;
Bode_O.title.String=' ';
Bode_O.Grid='on';
Bode_O.XLim={[1 1e4]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='off';

%% power stage parameter initialization
Tstart=0.015;
Rs=1e-1;
Vdcref=270;
alpha=1e10;         %% 0.5
Cdc=105e-6; 
R=96;
RCdc=0.049;
P=Vdcref^2/R;
Vse=57.5;
Vse=Vse+P/3/Vse*Rs*1;
Vsm=Vse*sqrt(2);
Vsm_o=Vse*sqrt(2);
L=470e-6;
RL=110e-3;
I = [1 0; 0 1];
f=60;%400;
omega=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;
Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd_s=Vsdq(1);
Vsq_s=Vsdq(2);
Id_ref=7.6307;
Iq_ref=0;
s=tf([1 0],[0 1]);

%% controller parameters calculation (need three sets of parameters)
kpi=Kpi_AFE/270;
kii=Kii_AFE/270;
kpv=Kpv_AFE;
kiv=Kiv_AFE;
tf_pll_z=tf([FWPI_a1 FWPI_a2],[1 -1],Tsw);
tf_pll=d2c(tf_pll_z);


%% signal conditioning filter
C1=75e-12;C2=39e-12;R1=15e3;R2=15e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));
tf_filter_dq = [tf_filter_dq(1,1) 0*tf_filter_dq(1,2); 0*tf_filter_dq(2,1) tf_filter_dq(2,2)];

Gdeldd = exp((-1/fsw)*s)*cos(-omega*(1.5/fsw));
Gdeldq = exp((-1/fsw)*s)*sin(-omega*(1.5/fsw));
Tdelay = 1.5/fsw;
Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);%exp(-0.5/fsw*s);
Gsvmdq = -(1-.5*Tdelay*s)*1/(1+.5*Tdelay*s);
Gsvm = [Gsvmdd 0*Gsvmdq;-0*Gsvmdq Gsvmdd];
Gsvm = tf(JF_DQFromABC(Gsvmdd,omega));
I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
Gdel = Gsvm;
Ki = tf_filter_dq;
Kv = tf_filter_dq;
fprintf('initialization is done!\n')

%% get the current loop gain:
%% Tid model linearization:
model = 'Ti_AFE';
io=getlinio(model);
%% Linearize the model
Tid = -linearize(model,0.3,io);
figHandle=figure(1);
set(figHandle, 'Position', [50, 10, 1000, 800]);
bode(Tid,Bode_O)
legend('Tid\_avg\_sim')
Bode_Darklines(3)
hold on
%% Tvdc model linearization:
model = 'Tv_AFE';
io=getlinio(model);
%% Linearize the model
Tvdc = -linearize(model,0.6,io);
figHandle=figure(2);
set(figHandle,'Position',[50, 10, 1000, 800]);
bode(Tvdc,Bode_O)
legend('Tvdc\_avg\_sim')
Bode_Darklines(3)