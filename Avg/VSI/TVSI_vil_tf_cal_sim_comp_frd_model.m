clc
clear all;
close all
%% voltage loop control parameters for different cases
Kpv_VSI      = 1*0.0495;%0.0158;%          % VSI voltage loop p controller           
Kiv_VSI      = 4*130.8997;%41.9;%        % VSI voltage loop i controller
Kpi_VSI      = 0.1*0.0465*270;      % VSI current loop p controller
Kii_VSI      = 0.0*121.0095*270;    % VSI current loop i controller
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
Bode_O.PhaseWrapping='off';

% bdclose all
% set_param(0,'CharacterEncoding','ISO-8859-1')
%% for phase to netural derivation
s=tf([1,0],[0,1]);
I = [tf([0 1],[0 1]) tf([0 0],[0 1]); tf([0 0],[0 1]) tf([0 1],[0 1])];
L=1.0e-3;%0.5e-3;%
RL=0.13;
R=12;
C=1*31.5e-6;
CR=.2;
Cdc=150e-6;
f=60;%400;%
omega=2*pi*f;
w=2*pi*f;
fsw=20e+3; Tsw=1/fsw; wsw=2*pi*fsw;

%% PWM delay
Gdeldd = exp((-1/fsw)*s)*cos(omega*(1/fsw));
Gdeldq = exp((-1/fsw)*s)*sin(omega*(1/fsw));
Tdelay = 1.5/fsw;%1.5/fsw;
Gsvmdd = (1-.5*Tdelay*s)/(1+.5*Tdelay*s);%exp(-0.5/fsw*s);
Gsvmdq = -(1-.5*Tdelay*s)*0/(1+.5*Tdelay*s);
Gsvm = tf(JF_DQFromABC(Gsvmdd,w));
Gsvm = [Gsvm(1,1) 0*Gsvm(1,2); -0*Gsvm(1,2) Gsvm(2,2)];
Gdel = Gsvm;
%% signal conditioning filter
C1=75e-12;C2=39e-12;R1=18e3;R2=18e3;
tf_filter = tf([0 1],[C1*C2*R1*R2 C2*(R1+R2) 1]);
tf_filter_dq = tf(JF_DQFromABC(tf_filter,omega));
tf_filter_dq = [tf_filter_dq(1,1) 0*tf_filter_dq(1,2); -0*tf_filter_dq(1,2) tf_filter_dq(2,2)];

Ki = tf_filter_dq;
Kv = tf_filter_dq;%Hv;

Vdc=270; Vdcref=270;
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsdq=[sqrt(3/2)*Vsm; 0];
Vdref=Vsdq(1); Vqref=Vsdq(2);

fv=1000;
fi=1000;
Kpi_VSI=1*2*pi*fi*L*0.5;
Kii_VSI=0*2*pi*fi*(RL)*10;
Kpv_VSI=2*pi*fv*C;
Kiv_VSI=2*pi*fv/R;


omega = f*2*pi;
Rac =  R;
Vod = Vdref;
Voq = Vqref;
Rl=RL;
Dd = 0.29997;
Dq = 0.08;
Idref = 8.2994;
Iqref = 7.8846;

kpi=Kpi_VSI/270;
kii=Kii_VSI/270;
kpv=Kpv_VSI;
kiv=Kiv_VSI;

%% get the current loop gain:
%% Tid model linearization:
model = 'Ti_VSI'; 
io=getlinio(model);
%% Linearize the model
Tid = -linearize(model,0.3,io);
figHandle=figure(1);
set(figHandle, 'Position', [50, 10, 1000, 800]);
bode(Tid,Bode_O)
legend('Tid\_avg\_sim')
Bode_Darklines(3)

%% Tvd model linearization:
model = 'Tv_VSI';
io=getlinio(model);
%% Linearize the model
Tvd = -linearize(model,0.6,io);
figHandle=figure(2);
set(figHandle,'Position',[50, 10, 1000, 800]);
bode(Tvd,Bode_O)
legend('Tvd\_avg\_sim')
Bode_Darklines(3)

