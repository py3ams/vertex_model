function thesis_graph_fig()

disp('busy');close all;

simulation_name = 'simulation_of_all_forces_with';
plot_name = 'angle_deviations';
save_plot_logical = 1;

linewidth = 2;

load(['Saves/',simulation_name,'/final_save.mat'])

figure('position',[100 100 325 300],'PaperPositionMode','auto','color','white')
set(gcf,'DefaultLineLineWidth',linewidth)

axes('position',[0.23 0.15 0.72 0.8])

statistics_counter = stats.counter;
time_range = linspace(0,1,stats.counter);

eval([plot_name,'_plot(linewidth,stats,statistics_counter,time_range);'])

% cell_area_plot(linewidth,stats,statistics_counter,time_range);

addpath('~/Documents/export_fig/')
if save_plot_logical
   if ~exist(['Figs/',simulation_name],'dir')
      mkdir('Figs/',simulation_name);
   end
%    saveas(gcf,['Figs/',simulation_name,'/',simulation_name,'_',plot_name,'.eps'],'psc2')
   export_fig(['Figs/',simulation_name,'/',simulation_name,'_',plot_name,'.eps'],'-nocrop');
end

end

function height_to_area_ratios_plot(linewidth,stats,statistics_counter,time_range)

cell_height_to_area = stats.cell_height_to_area(1:statistics_counter,:);
plot(time_range,cell_height_to_area(:,1))
hold on
plot(time_range,cell_height_to_area(:,1)+2*cell_height_to_area(:,4),'r')
% plot(time_range,cell_height_to_area(:,1)-2*cell_height_to_area(:,4),'r')

set(gca,'FontName','arial','fontweight','bold','fontsize',13);
xlabel('Time')
ylabel('Height area ratio')

axis([0 1 0 35])

set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))

end

function angle_deviations_plot(linewidth,stats,statistics_counter,time_range)

angle_deviation = stats.angle_deviation(1:statistics_counter,:);
plot(time_range,angle_deviation(:,1))
set(gca,'FontName','arial','fontweight','bold','fontsize',13);
xlabel('Time')
ylabel('Angle deviation')

axis([0 1 0 0.4])

set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
set(gca,'YTickLabel',sprintf('%0.1f|',str2num(get(gca,'YTickLabel'))))

end

function edge_lengths_plot(linewidth,stats,statistics_counter,time_range)

edge_length = stats.edge_length(1:statistics_counter,:);
size(stats.edge_length)
plot(time_range,edge_length(:,1))
hold on
% plot(time_range,edge_length(:,2),'r')
% plot(time_range,edge_length(:,3),'r')
plot(time_range,edge_length(:,1)+2*edge_length(:,4),'r')
plot(time_range,edge_length(:,1)-2*edge_length(:,4),'r')

set(gca,'FontName','arial','fontweight','bold','fontsize',13);
xlabel('Time')
ylabel('Edge length')

axis([0 1 0 0.15])

set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))

end

function plot_name = cell_perimeter_plot(linewidth,stats,statistics_counter,time_range)

plot_name = 'cell_perimeters';

cell_perimeter = stats.cell_perimeter(1:statistics_counter,:);
plot(time_range,cell_perimeter(:,1))
hold on
plot(time_range,cell_perimeter(:,1)+2*cell_perimeter(:,4),'r')
plot(time_range,cell_perimeter(:,1)-2*cell_perimeter(:,4),'r')
set(gca,'FontName','arial','fontweight','bold','fontsize',13);
xlabel('Time')
ylabel('Cell perimeter')

axis([0 1 0.3 0.451])

set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))

end

function plot_name = boundary_force_plot(linewidth,stats,statistics_counter,time_range)

plot_name = 'boundary_forces_per_vertex';

total_boundary_deformation_force =...
	stats.total_boundary_deformation_force(1:statistics_counter,:);
total_boundary_edge_force = stats.total_boundary_edge_force(1:statistics_counter,:);
total_boundary_vertices = stats.total_boundary_vertices(1:statistics_counter,:);

plot(time_range,total_boundary_deformation_force./total_boundary_vertices)
hold all
plot(time_range,total_boundary_edge_force./total_boundary_vertices)
set(gca,'FontName','arial','fontweight','bold','fontsize',13);
xlabel('Time')
ylabel('Mean force per vertex')
legend('C_{B1}','C_{B2}','location','best')

axis([0 1 0 0.151])

set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))

end

function plot_name = force_plot(linewidth,stats,statistics_counter,time_range)

plot_name = 'forces_per_vertex';

total_area_force = stats.total_area_force(1:statistics_counter,:);
total_deformation_force = stats.total_deformation_force(1:statistics_counter,:);
total_elongation_force = stats.total_elongation_force(1:statistics_counter,:);
total_perimeter_force = stats.total_perimeter_force(1:statistics_counter,:);
total_tension_force = stats.total_tension_force(1:statistics_counter,:);
total_no_vertices = stats.total_no_vertices(1:statistics_counter,:);
plot(time_range,total_area_force./total_no_vertices)
hold all
plot(time_range,total_deformation_force./total_no_vertices)
plot(time_range,total_elongation_force./total_no_vertices)
plot(time_range,total_perimeter_force./total_no_vertices)
plot(time_range,total_tension_force./total_no_vertices)
set(gca,'FontName','arial','fontweight','bold','fontsize',13);
xlabel('Time')
ylabel('Mean force per vertex')
legend('C_A','C_D','C_H','C_P','C_L','location','best')

axis([0 1 0 0.04])

set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))

end

function cell_areas_plot(linewidth,stats,statistics_counter,time_range)

cell_area = stats.cell_area(1:statistics_counter,:);

plot(time_range,cell_area(:,1),'linewidth',linewidth)
hold on
plot(time_range,cell_area(:,1)+2*cell_area(:,4),'r','linewidth',linewidth)
plot(time_range,cell_area(:,1)-2*cell_area(:,4),'r','linewidth',linewidth)
set(gca,'FontName','arial','fontweight','bold','fontsize',13);
xlabel('Time')
ylabel('Cell area')

axis([0 1 0.002 0.016])

set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
set(gca,'YTickLabel',sprintf('%0.3f|',str2num(get(gca,'YTickLabel'))*0.001))

end



