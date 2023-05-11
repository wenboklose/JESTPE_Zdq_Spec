close all
clear all
%% Bode plot options
Bode_O=bodeoptions;
Bode_O.XLabel.FontSize=14;
Bode_O.YLabel.FontSize=14;
Bode_O.TickLabel.FontSize=14;
Bode_O.Title.FontSize=14;
Bode_O.title.String=' ';
Bode_O.Grid='on';
Bode_O.XLim={[1 1e3]};
Bode_O.XLimMode={'manual'};
Bode_O.FreqUnits='Hz';
Bode_O.PhaseWrapping='on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linewidth=3; fontsize=18;
linestyle1='b'; linestyle2='g';
linestyle3='r'; linestyle4='m';
linestyle5='k'; linestyle6='c';

linestyle=[linestyle1; linestyle2; linestyle3; linestyle4];
s=tf([1 0],[0 1]);
load ZVSI_v1.mat
Z_S = ZVSI_ac_i_v_dec;
Z_S_dd = Z_S(1,1);
Z_S_qq = Z_S(2,2);
load ZAFE_96.mat
Z_AFE = Zac_AFE_i_dec;
C=60e-6;
YC=[C*s -C*60*2*pi;C*60*2*pi C*s];
L=200e-3;
ZL=[L*s -L*60*2*pi;L*60*2*pi L*s];
Z_L = 1/(1/Z_AFE+0*YC);
Z_L_dd = Z_L(1,1);
Z_L_qq = Z_L(2,2);
LL=(Z_S/Z_L);
Ldd = (Z_S_dd/Z_L_dd);
Lqq = (Z_S_qq/Z_L_qq);


Z_S_s = [Z_S_dd 0; 0 Z_S_qq];
Y_L_s = [1/Z_L_dd 0; 0 1/Z_L_qq];

% RR is the full GNC return ration
RR=Z_S/Z_L;
% rr is the return ratio defined in Rolando's paper
rr=Z_S_s*Y_L_s;
ll=Z_S*Y_L_s;
RRdd=ll(1,1);
RRdq=ll(1,2);
RRqd=ll(2,1);
RRqq=ll(2,2);
% calculate ACindex defined by Rolando
ACindex=(RRdd-RRqq)*(RRdd-RRqq)/(4*RRqd*RRdq);

% Eigenvalues and Sorting

% plot GNC for using RR and L
    n=1e3;
    f=logspace(-1,2.5,n);
%     f=1:1:1e4;
    w=f*2*pi;
    % For negative frequencies turn on
%     if(1)
%         w=[-fliplr(w) w];
%     %     n=n*2;
%     end
    % L(s) frequency response
    Lresp=freqresp(RR,w);
    lresp=freqresp(rr,w)
    for k=1:length(w)
%         Lddeigenvalues(:,k)=eig(Lddresp(:,:,k));
%         Lqqeigenvalues(:,k)=eig(Lqqresp(:,:,k));
        Leigenvalues(:,k)=eig(Lresp(:,:,k));
        leigenvalues(:,k)=eig(lresp(:,:,k));
%         leigenvalues(:,k)=eig(rr_res(:,:,k));
    end
%     [Leigenvalues]=sortloci(Leigenvalues);
%     figure(7);
%     plot(1./Leigenvalues(1,:),linestyle(1),'LineWidth',linewidth)
%     hold on
%     plot(1./Leigenvalues(2,:),strcat(':',linestyle(2)),'LineWidth',linewidth)
    
    
    fig=figure(7);
    set(fig, 'Position', [50, 10, 1000, 800]);
    plot(Leigenvalues(1,:),linestyle(1),'LineWidth',linewidth)
    hold on
    plot(Leigenvalues(2,:),linestyle(2),'LineWidth',linewidth)
    hold on
    plot(leigenvalues(1,:),strcat(':',linestyle(3)),'LineWidth',linewidth)
    hold on
    plot(leigenvalues(2,:),strcat(':',linestyle(4)),'LineWidth',linewidth)
    grid on
    axisloci=axis;
    legend({'{\it\lambda}_{1}','{\it\lambda}_{2}','Zsdd/Zldd','Zsqq/Zlqq'},'Fontsize',fontsize,'FontWeight','bold')
    plot(-1,0,'r+','LineWidth',linewidth)
%     hold off
    title('Characteristic Loci of L','Fontsize',fontsize,'FontWeight','bold')
    grid on
    set(gca,'FontSize',fontsize);
    
    h=ezplot('x^2+y^2=1');
    set(h,'color','r');
    
    figure(8)
    pzmap(LL)
    hold on
    
    figure(9)
    bode(Z_S,Z_L,Bode_O)
    Bode_Darklines(3)
    hold on
    Z_S_dd=Z_S(1,1);
    Z_L_dd=Z_L(1,1);
    Z_S_qq=Z_S(2,2);
    Z_L_qq=Z_L(2,2);
    figure(10)
    bode(Z_S_qq,Z_L_qq,Bode_O)
    Bode_Darklines(3)
    hold on
    figure(11)
    nyquist(Z_S_qq/Z_L_qq)
%     Bode_Darklines(3)
    hold on
    figure(12)
    bode(Z_S_dd,Z_L_dd)
    Bode_Darklines(3)
    hold on
    figure(13)
    nyquist(Z_S_dd/Z_L_dd)
    hold on
    
    figure(14)
    nyquist(Z_S)
    figure(15)
    nyquist(1/Z_L)
    Bode_Darklines(3)