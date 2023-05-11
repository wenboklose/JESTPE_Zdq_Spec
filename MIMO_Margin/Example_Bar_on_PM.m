%% verification of Bar-on example:
clc
clear all
close all
k_d = -7;
k_g = 1;

A=[0 1 0 0;-3 -0.75 1 0.25;0 0 0 1;4 1 -4 -1];
B=[0 0;0 1;0 0;0.25 0];
C=[1 0 0 0;0 0 1 0];
D=0;
sys=tf(ss(A,B,C,D));
s=tf('s');

% system 1
den = s^4+1.75*s^3+7.5*s^2+4*s+8;
num11 = 0.0625*s+0.25;
num12 = s^2+s+4;

num21 = 0.25*s^2+0.1875*s+0.75;
num22 = s+4;
L = sys;
L = 1/den*[num11 num12; num21 num22];
I = [tf([0 1],[0 1]), 0; 0 tf([0 1],[0 1])];
PM=-60;
U=[exp(-1i*PM*pi/180*s) 0;0 exp(-1i*PM*pi/180*s)];
M=L*U/(I+L*U);
figure
step(M)
% figure(1)
% nyquist(L)
% grid on
% hold on
% nyquist(M)
% freq_r = logspace(-5, 5,5000)*2*pi;
% Resp_L = freqresp(L,freq_r);
% Resp_M = freqresp(M,freq_r);
% Resp_I = freqresp(I,freq_r);
% for k=1:length(freq_r)
%    eig_L(:,k) = eig(Resp_L(:,:,k));             % eigenvalue of L
%    eig_L_I(:,k) = eig(Resp_L(:,:,k));           % calculate eigenvalue of L
%    [V,D]=eig(Resp_L(:,:,k));                    % calculate eigenvalue and eigenvector of L
%    C_V(k)=dot(V(:,1),V(:,2));                  % dot product of eigenvectors, test whether L is normal matrix
%    [V,D]=eig(Resp_L(:,:,k)'*Resp_L(:,:,k));
%    C_V_s(k)=dot(V(:,1),V(:,2));
%    Singular_L(:,k) = svd(Resp_L(:,:,k));        % calculate singular value of L
%    [U H]=poldec(Resp_L(:,:,k));                 % polar decomposatoin of L
%    eig_U(:,k)=eig(U);                           % calculate eigenvalue of Unitary part of L
%    eig_H(:,k)=eig(H);                           % calculate eigenvalue of Hermitain part of L
%    Singular_M(:,k) = svd(Resp_M(:,:,k));
%    [U H]=poldec(Resp_M(:,:,k));
%    eig_U_M(:,k)=eig(U);
%    eig_H_M(:,k)=eig(H);
% end
% figure(2) %% show the singular value of L is the same of the eigenvalue of L's Hermitain part, they are the boundary of L's eigenvalue
% semilogx(freq_r,Singular_L(1,:),'b*');
% hold on
% semilogx(freq_r,Singular_L(2,:),'b*');
% semilogx(freq_r,eig_H(1,:),'r-');
% semilogx(freq_r,eig_H(2,:),'r-');
% semilogx(freq_r,abs(eig_L(1,:)),'k');
% semilogx(freq_r,abs(eig_L(2,:)),'k');
% grid on
% figure(3) %% show the same thing for M as been shown for L in figure(2)
% semilogx(freq_r,Singular_M(1,:),'b*');
% hold on
% semilogx(freq_r,Singular_M(2,:),'b*');
% semilogx(freq_r,eig_H_M(1,:),'r-');
% semilogx(freq_r,eig_H_M(2,:),'r-');
% 
% semilogx(freq_r,Singular_L(1,:),'g*');
% semilogx(freq_r,Singular_L(2,:),'g*');
% 
% 
% grid on
% figure(4) %% show the eigenvalue of L in complex plane
% plot(eig_L(1,1:length(eig_L)),'b')
% hold on
% plot(eig_L(2,1:length(eig_L)),'r')
% grid on
% ezplot('x^2+y^2=1')
% Mag=abs(eig_L(2,:));
% Phase=angle(eig_L(2,:))*180/pi;
% figure(41) %% show magnitude and phase of L's eigenvalue (first one)
% subplot(2,1,1)
% semilogx(freq_r,20*log10(Mag),'b*');
% grid on
% subplot(2,1,2)
% semilogx(freq_r,Phase,'r*');
% hold off
% grid on
% Mag=abs(eig_L(1,:));
% Phase=angle(eig_L(1,:))*180/pi;
% figure(411) %% show magnitude and phase of L's eigenvalue (second one)
% subplot(2,1,1)
% semilogx(freq_r,20*log10(Mag),'b*');
% grid on
% subplot(2,1,2)
% semilogx(freq_r,Phase,'r*');
% hold off
% grid on
% 
% %% find the biggest singular value of M
% [M_max_1 index_1]=max(Singular_M(1,:));
% [M_max_2 index_2]=max(Singular_M(2,:));
% if M_max_1>=M_max_2
%     M_max = M_max_1;
%     GM = 1/M_max_1;
%     freq_crit = freq_r(index_1); 
%     index = index_1;
% else
%     M_max = M_max_2;
%     GM = 1/M_max_2;
%     freq_crit = freq_r(index_2);
%     index = index_2;
% end
% GM
% Resp_M(:,:,index);
% [U,S,V]=svd(Resp_M(:,:,index));
% % [U,S,V]=svd([0.0508-0.0777j -0.1881-1.8974j;0.0921-0.2726j -2.0354-3.1041j])
% um1=U(:,1);vm1=V(:,1);
% PM=acosd(real(-um1'*vm1))
% 
% 
% % GM = 1/M_max_2
% %% now, form a perturbation which will destabilize the system
% D = -1*GM*vm1*um1';
% % D=[0.0822+0.1281*j -0.0125-0.0960*j;0.0141-0.1191*j 0.1122-0.1733*j];
% [U R]=poldec(D);
% eig(R);
% [S V]=eig(U);
% angle(eig(U))*180/pi;
% max(abs(angle(eig(U))))*180/pi;
% 
% figure(42)
% plot(eig_L_I(1,1:length(eig_L)),'b')
% hold on
% plot(eig_L_I(2,1:length(eig_L)),'r')
% grid on
% plot(conj(eig_L_I(1,1:length(eig_L))),'b+')
% plot(conj(eig_L_I(2,1:length(eig_L))),'r+')
% for k=1:length(freq_r)
%    eig_L(:,k) = eig(Resp_L(:,:,k)*(Resp_I(:,:,k)+D));
% end
% [eig_L]=sortloci(eig_L);
% 
% plot(eig_L(1,1:length(eig_L)),'r')
% hold on
% plot(eig_L(2,1:length(eig_L)),'b')
% grid on
% plot(conj(eig_L(1,1:length(eig_L))),'r*')
% plot(conj(eig_L(2,1:length(eig_L))),'b*')
% 
% %% this part is for multiply perturbation: trying to define phase and gain margin using eigenvalue of matrix L
% % first, get the magnitude and angle of L's eigenvalue:
% for k=1:length(freq_r)
%    eig_L(:,k) = eig(Resp_L(:,:,k));
% end
% 
% angle_eigs_L=angle(eig_L);
% Mag_eigs_L=abs(eig_L);
% 
% % for phase margin, find at which frequency Mag of L' eig reach unity
% % circle
% 
% index_1=find(Mag_eigs_L(1,:)>0.999&Mag_eigs_L(1,:)<1.002);
% index_2=find(Mag_eigs_L(2,:)>0.99&Mag_eigs_L(2,:)<1.01);
% PM_1=(angle_eigs_L(1,index_1)+pi)*180/pi;
% PM_2=(angle_eigs_L(2,index_2)+pi)*180/pi;
% 
% if isempty(PM_1)
%     [PM Index]= min(PM_2);
%     Index = index_2;
% elseif isempty(PM_2)
%     [PM Index] = min(PM_1);
%     Index = index_1;
% elseif min(PM_1)>=min(PM_2)
%     [PM Index] = min(PM_2);
%     Index = index_2;
% else
%     [PM Index] = min(PM_1);
%     Index = index_1;
% end
% Index
% freq_r(Index)/(2*pi);
% Phase_Margin=PM
% eig_L(:,Index);
% for k=1:length(freq_r)
%    eig_L(:,k) = eig(Resp_L(:,:,k));
% end
% % [eig_L]=sortloci(eig_L);
% figure(5)
% plot(eig_L(1,1:length(eig_L)),'b')
% hold on
% plot(eig_L(2,1:length(eig_L)),'r')
% grid on
% ezplot('x^2+y^2=1')
% 
% U=[exp(-1i*PM*pi/180) 0;0 exp(-1i*PM*pi/180)];
% for k=1:length(freq_r)
%    eig_L_1(:,k) = eig(Resp_L(:,:,k)*(U));
% end
% [eig_L_1]=sortloci(eig_L_1);
% % figure(41)
% plot(eig_L_1(1,1:length(eig_L_1)),'b*')
% % hold on
% plot(eig_L_1(2,1:length(eig_L_1)),'r*')
% % grid on
% 
% U_1=[exp(-1i*(PM-k_d)*pi/180) 0;0 exp(-1i*(PM-k_d)*pi/180)];
% figure(51)
% for k=1:length(freq_r)
%    eig_L_2(:,k) = eig(Resp_L(:,:,k)*U_1);
% end
% [eig_L_2]=sortloci(eig_L_2);
% U_2=[exp(-1i*(PM+k_d)*pi/180) 0;0 exp(-1i*(PM+k_d)*pi/180)];
% figure(51)
% for k=1:length(freq_r)
%    eig_L_3(:,k) = eig(Resp_L(:,:,k)*U_2);
% end
% [eig_L_3]=sortloci(eig_L_3);
% % figure(41)
% plot(eig_L(1,1:length(eig_L)),'b')
% hold on
% plot(eig_L(2,1:length(eig_L)),'r')
% 
% plot(eig_L_1(1,1:length(eig_L_1)),'b*')
% hold on
% plot(eig_L_1(2,1:length(eig_L_1)),'r*')
% grid on
% 
% plot(eig_L_2(1,1:length(eig_L_2)),'bo')
% hold on
% plot(eig_L_2(2,1:length(eig_L_2)),'ro')
% grid on
% 
% plot(eig_L_3(1,1:length(eig_L_2)),'bo')
% hold on
% plot(eig_L_3(2,1:length(eig_L_2)),'ro')
% grid on
% ezplot('x^2+y^2=1')
% % for gain margin, find at which frequency phase of L' eig reach -pi
% 
% index_1=find(angle_eigs_L(1,:)>-pi-0.02&angle_eigs_L(1,:)<-pi+0.02);
% index_2=find(angle_eigs_L(2,:)>-pi-0.02&angle_eigs_L(2,:)<-pi+0.02);
% GM_1=(1./Mag_eigs_L(1,index_1));
% GM_2=(1./Mag_eigs_L(2,index_2));
% 
% if isempty(GM_1)
%     [GM Index] = min(GM_2);
%     Index = index_2;
% elseif isempty(GM_2)
%     [GM Index] = min(GM_1);
%     Index = index_1;
% elseif min(GM_1)>=min(GM_2)
%     [GM Index] = min(GM_2);
%     Index = index_2;
% else
%     [GM Index] = min(GM_1);
%     Index = index_1;
% end
% Index;
% freq_r(Index)/(2*pi);
% Gain_Margin=20*log10(GM)
% eig_L(:,Index);
% for k=1:length(freq_r)
%    eig_L(:,k) = eig(Resp_L(:,:,k));
% end
% % [eig_L]=sortloci(eig_L);
% figure(6)
% plot(eig_L(1,1:length(eig_L)),'b')
% hold on
% plot(eig_L(2,1:length(eig_L)),'r')
% grid on
% ezplot('x^2+y^2=1')
% R=[GM 0;0 GM];
% for k=1:length(freq_r)
%    eig_L_1(:,k) = eig(Resp_L(:,:,k)*(R));
% end
% [eig_L_1]=sortloci(eig_L_1);
% % figure(41)
% plot(eig_L_1(1,1:length(eig_L_1)),'b*')
% hold on
% plot(eig_L_1(2,1:length(eig_L_1)),'r*')
% grid on
% 
% R_1=[k_g*GM 0;0 k_g*GM];
% % step response:
% L_2=L*R_1;
% M_2=L_2/(I+L_2);
% figure(66)
% step(M)
% hold on
% step(M_2)
% grid on
% 
% R_2=[1/k_g*GM 0;0 1/k_g*GM];
% % step response:
% L_3=L*R_2;
% M_3=L_3/(I+L_3);
% figure(661)
% step(M)
% hold on
% step(M_3)
% grid on
% 
% 
% figure(61)
% for k=1:length(freq_r)
%    eig_L_2(:,k) = eig(Resp_L(:,:,k)*R_1);
% end
% [eig_L_2]=sortloci(eig_L_2);
% 
% R_2=[1/k_g*GM 0;0 1/k_g*GM];
% for k=1:length(freq_r)
%    eig_L_3(:,k) = eig(Resp_L(:,:,k)*R_2);
% end
% [eig_L_2]=sortloci(eig_L_2);
% 
% % figure(41)
% plot(eig_L(1,1:length(eig_L)),'b')
% hold on
% plot(eig_L(2,1:length(eig_L)),'r')
% 
% plot(eig_L_1(1,1:length(eig_L_1)),'b*')
% hold on
% plot(eig_L_1(2,1:length(eig_L_1)),'r*')
% grid on
% 
% plot(eig_L_2(1,1:length(eig_L_2)),'bo')
% hold on
% plot(eig_L_2(2,1:length(eig_L_2)),'ro')
% grid on
% 
% plot(eig_L_3(1,1:length(eig_L_2)),'bo')
% hold on
% plot(eig_L_3(2,1:length(eig_L_2)),'ro')
% grid on
% ezplot('x^2+y^2=1')
% 
% 
% 
% 
% %% this part is for multiply perturbation: trying to define phase and gain margin using singular value of matrix L
% % first, get the magnitude and angle of L's eigenvalue:
% for k=1:length(freq_r)
%    eig_L(:,k) = eig(Resp_L(:,:,k));
% end
% 
% angle_eigs_L=angle(eig_U);
% Mag_eigs_L=abs(eig_H);
% 
% % for phase margin, find at which frequency Mag of L' eig reach unity
% % circle
% 
% index_1=find(Mag_eigs_L(1,:)>0.999&Mag_eigs_L(1,:)<1.002);
% index_2=find(Mag_eigs_L(2,:)>0.99&Mag_eigs_L(2,:)<1.01);
% PM_1=(angle_eigs_L(1,index_1)+pi)*180/pi;
% PM_2=(angle_eigs_L(2,index_2)+pi)*180/pi;
% 
% if isempty(PM_1)
%     [PM Index]= min(PM_2);
%     Index = index_2;
% elseif isempty(PM_2)
%     [PM Index] = min(PM_1);
%     Index = index_1;
% elseif min(PM_1)>=min(PM_2)
%     [PM Index] = min(PM_2);
%     Index = index_2;
% else
%     [PM Index] = min(PM_1);
%     Index = index_1;
% end
% Index
% freq_r(Index)/(2*pi);
% Phase_Margin=PM
% eig_L(:,Index);
% for k=1:length(freq_r)
%    eig_L(:,k) = eig(Resp_L(:,:,k));
% end
% % [eig_L]=sortloci(eig_L);
% figure(7)
% plot(eig_L(1,1:length(eig_L)),'b')
% hold on
% plot(eig_L(2,1:length(eig_L)),'r')
% grid on
% ezplot('x^2+y^2=1')
% 
% U=[exp(-1i*PM*pi/180) 0;0 exp(-1i*PM*pi/180)];
% for k=1:length(freq_r)
%    eig_L_1(:,k) = eig(Resp_L(:,:,k)*(U));
% end
% [eig_L_1]=sortloci(eig_L_1);
% % figure(41)
% plot(eig_L_1(1,1:length(eig_L_1)),'b*')
% % hold on
% plot(eig_L_1(2,1:length(eig_L_1)),'r*')
% 
% 
% 
% 
% % grid on
% U_1=[exp(-1*(PM-k_d)*pi/180*s/freq_r(index)) 0;0 exp(-1*(PM-k_d)*pi/180*s/freq_r(index))];
% % step response:
% L_1=L*U_1;
% M_1=L_1/(I+L_1);
% figure(77)
% step(M)
% hold on
% step(M_1)
% grid on
% 
% 
% U_2=[exp(-1*(PM+k_d)*pi/180*s/freq_r(index)) 0;0 exp(-1*(PM+k_d)*pi/180*s/freq_r(index))];
% 
% 
% % step response:
% L_2=L*U_2;
% M_2=L_2/(I+L_2);
% figure(771)
% step(M)
% hold on
% step(M_2)
% grid on
% 
% 
% U_1=[exp(-1i*(PM-k_d)*pi/180) 0;0 exp(-1i*(PM-k_d)*pi/180)];
% figure(71)
% for k=1:length(freq_r)
%    eig_L_2(:,k) = eig(Resp_L(:,:,k)*U_1);
% end
% [eig_L_2]=sortloci(eig_L_2);
% U_2=[exp(-1i*(PM+k_d)*pi/180) 0;0 exp(-1i*(PM+k_d)*pi/180)];
% figure(71)
% for k=1:length(freq_r)
%    eig_L_3(:,k) = eig(Resp_L(:,:,k)*U_2);
% end
% [eig_L_3]=sortloci(eig_L_3);
% % figure(41)
% plot(eig_L(1,1:length(eig_L)),'b')
% hold on
% plot(eig_L(2,1:length(eig_L)),'r')
% 
% plot(eig_L_1(1,1:length(eig_L_1)),'b*')
% hold on
% plot(eig_L_1(2,1:length(eig_L_1)),'r*')
% grid on
% 
% plot(eig_L_2(1,1:length(eig_L_2)),'bo')
% hold on
% plot(eig_L_2(2,1:length(eig_L_2)),'ro')
% grid on
% 
% plot(eig_L_3(1,1:length(eig_L_2)),'bo')
% hold on
% plot(eig_L_3(2,1:length(eig_L_2)),'ro')
% grid on
% ezplot('x^2+y^2=1')
% % for gain margin, find at which frequency phase of L' eig reach -pi
% 
% % index_1=find(angle_eigs_L(1,:)>-pi-0.02&angle_eigs_L(1,:)<-pi+0.02);
% % index_2=find(angle_eigs_L(2,:)>-pi-0.02&angle_eigs_L(2,:)<-pi+0.02);
% % GM_1=(1./Mag_eigs_L(1,index_1));
% % GM_2=(1./Mag_eigs_L(2,index_2));
% % 
% % if isempty(GM_1)
% %     [GM Index] = min(GM_2);
% %     Index = index_2;
% % elseif isempty(GM_2)
% %     [GM Index] = min(GM_1);
% %     Index = index_1;
% % elseif min(GM_1)>=min(GM_2)
% %     [GM Index] = min(GM_2);
% %     Index = index_2;
% % else
% %     [GM Index] = min(GM_1);
% %     Index = index_1;
% % end
% % Index
% % freq_r(Index)/(2*pi)
% % Gain_Margin=20*log10(GM)
% % eig_L(:,Index);
% % for k=1:length(freq_r)
% %    eig_L(:,k) = eig(Resp_L(:,:,k));
% % end
% % % [eig_L]=sortloci(eig_L);
% % figure(8)
% % plot(eig_L(1,1:length(eig_L)),'b')
% % hold on
% % plot(eig_L(2,1:length(eig_L)),'r')
% % grid on
% % ezplot('x^2+y^2=1')
% % R=[GM 0;0 GM];
% % for k=1:length(freq_r)
% %    eig_L_1(:,k) = eig(Resp_L(:,:,k)*(R));
% % end
% % [eig_L_1]=sortloci(eig_L_1);
% % % figure(41)
% % plot(eig_L_1(1,1:length(eig_L_1)),'b*')
% % hold on
% % plot(eig_L_1(2,1:length(eig_L_1)),'r*')
% % grid on
% % 
% % R_1=[k_g*GM 0;0 k_g*GM];
% % % step response:
% % L_2=L*R_1;
% % M_2=L_2/(I+L_2);
% % figure(88)
% % step(M)
% % hold on
% % step(M_2)
% % grid on
% % 
% % R_2=[1/k_g*GM 0;0 1/k_g*GM];
% % % step response:
% % L_3=L*R_2;
% % M_3=L_3/(I+L_3);
% % figure(881)
% % step(M)
% % hold on
% % step(M_3)
% % grid on
% % 
% % 
% % figure(81)
% % for k=1:length(freq_r)
% %    eig_L_2(:,k) = eig(Resp_L(:,:,k)*R_1);
% % end
% % [eig_L_2]=sortloci(eig_L_2);
% % 
% % R_2=[1/k_g*GM 0;0 1/k_g*GM];
% % for k=1:length(freq_r)
% %    eig_L_3(:,k) = eig(Resp_L(:,:,k)*R_2);
% % end
% % [eig_L_2]=sortloci(eig_L_2);
% % 
% % % figure(41)
% % plot(eig_L(1,1:length(eig_L)),'b')
% % hold on
% % plot(eig_L(2,1:length(eig_L)),'r')
% % 
% % plot(eig_L_1(1,1:length(eig_L_1)),'b*')
% % hold on
% % plot(eig_L_1(2,1:length(eig_L_1)),'r*')
% % grid on
% % 
% % plot(eig_L_2(1,1:length(eig_L_2)),'bo')
% % hold on
% % plot(eig_L_2(2,1:length(eig_L_2)),'ro')
% % grid on
% % 
% % plot(eig_L_3(1,1:length(eig_L_2)),'bo')
% % hold on
% % plot(eig_L_3(2,1:length(eig_L_2)),'ro')
% % grid on
% % ezplot('x^2+y^2=1')