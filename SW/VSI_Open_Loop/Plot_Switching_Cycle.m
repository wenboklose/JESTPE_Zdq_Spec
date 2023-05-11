figHandle=figure
set(figHandle,'Position',[50, 10, 1000, 800]);

subplot(6,1,1)
plot(gate_vab.time,gate_vab.signals(1).values)
xlim([101/(20e3) 102/(20e3)])
ylim([0 1.1])
% xlabel('Time (s)','FontSize',14,'FontWeight','normal','Color','k')
% ylabel('Current (s)','FontSize',14,'FontWeight','normal','Color','k')
title('Sa','FontSize',14)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',3)
set(gca,'FontSize',14)
grid on

subplot(6,1,2)
plot(gate_vab.time,gate_vab.signals(2).values)
xlim([101/(20e3) 102/(20e3)])
ylim([0 1.1])
% xlabel('Time (s)','FontSize',14,'FontWeight','normal','Color','k')
% ylabel('Current (s)','FontSize',14,'FontWeight','normal','Color','k')
title('Sb','FontSize',14)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',3)
set(gca,'FontSize',14)
grid on

subplot(6,1,3)
plot(gate_vab.time,gate_vab.signals(3).values)
xlim([101/(20e3) 102/(20e3)])
ylim([0 1.1])
% xlabel('Time (s)','FontSize',14,'FontWeight','normal','Color','k')
% ylabel('Current (s)','FontSize',14,'FontWeight','normal','Color','k')
title('Sc','FontSize',14)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',3)
set(gca,'FontSize',14)
grid on

subplot(6,1,4)
plot(gate_vab.time,gate_vab.signals(4).values)
xlim([101/(20e3) 102/(20e3)])
ylim([-300 300])
% xlabel('Time (s)','FontSize',14,'FontWeight','normal','Color','k')
ylabel('Voltage (V)','FontSize',14,'FontWeight','normal','Color','k')
title('Vab = (Sa-Sb)*Vdc','FontSize',14)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',3)
set(gca,'FontSize',14)
grid on

subplot(6,1,5)
plot(gate_vab.time,gate_vab.signals(5).values)
xlim([101/(20e3) 102/(20e3)])
ylim([-300 300])
% xlabel('Time (s)','FontSize',14,'FontWeight','normal','Color','k')
ylabel('Voltage (V)','FontSize',14,'FontWeight','normal','Color','k')
title('Vbc = (Sb-Sc)*Vdc','FontSize',14)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',3)
set(gca,'FontSize',14)
grid on

subplot(6,1,6)
plot(gate_vab.time,gate_vab.signals(6).values)
xlim([101/(20e3) 102/(20e3)])
ylim([-300 300])
xlabel('Time (s)','FontSize',14,'FontWeight','normal','Color','k')
ylabel('Voltage (V)','FontSize',14,'FontWeight','normal','Color','k')
title('Vca = (Sc-Sa)*Vdc','FontSize',14)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',3)
set(gca,'FontSize',14)
grid on