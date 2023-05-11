figHandle=figure
set(figHandle,'Position',[50, 10, 1000, 800]);
subplot(2,1,1)
plot(I_IGBT.time,I_IGBT.signals.values)
grid on
xlim([0.01-2/400 0.01])
xlabel('Time (s)','FontSize',18,'FontWeight','normal','Color','k')
ylabel('Current (A)','FontSize',18,'FontWeight','normal','Color','k')
title('IGBT Current (Phase A) 9.85 A RMS','FontSize',18)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',1)
set(gca,'FontSize',18)

subplot(2,1,2)
plot(Ia_L.time,Ia_L.signals.values)
xlim([0.01-2/400 0.01])
xlabel('Time (s)','FontSize',18,'FontWeight','normal','Color','k')
ylabel('Current (A)','FontSize',18,'FontWeight','normal','Color','k')
title('Inductor Current (Phase A) 13.85 A RMS','FontSize',18)
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',1)
set(gca,'FontSize',18)
grid on