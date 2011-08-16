disp('busy');clear all;close all;

simulation_name = 'deformation_force_only';
plot_name = 'angle_deviations';

save_plot_logical = 1;
linewidth = 2;

load(['Saves/',simulation_name,'/final_save.mat'])
time_range = linspace(0,1,stats.counter);

figure('position',[100 100 300 300],'PaperPositionMode','auto')

angle_deviation = stats.angle_deviation(1:stats.counter,:);
plot(time_range,angle_deviation(:,1),'linewidth',linewidth)
hold on
% plot(time_range,angle_deviation(:,2),'r')
% plot(time_range,angle_deviation(:,3),'r')
% plot(time_range,angle_deviation(:,1)+2*angle_devi5ation(:,4),'r','linewidth',linewidth)
% plot(time_range,angle_deviation(:,1)-2*angle_deviation(:,4),'r','linewidth',linewidth)
set(gca,'FontName','arial','fontweight','bold','fontsize',13);
xlabel('Time')
ylabel('Angle deviation')
% axis([0 time min(angle_deviation(:,1)-2.5*angle_deviation(:,4)) ...
% 		max(angle_deviation(:,1)+2.5*angle_deviation(:,4))])

axis([0 1 0 0.4])

set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
set(gca,'YTickLabel',sprintf('%0.1f|',str2num(get(gca,'YTickLabel'))))

if save_plot_logical
   if ~exist(['Figs/',simulation_name],'dir')
      mkdir('Figs/',simulation_name);
   end
   saveas(gcf,['Figs/',simulation_name,'/',simulation_name,'_',plot_name,'.eps'],'psc2')
end