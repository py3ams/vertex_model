function thesis_graph_fig()

disp('busy');close all;

folder_name = 'simulation_of_all_forces';
plot_name = 'boundary_forces';
save_plot_logical = 1;

linewidth = 2;

load(['Saves/',folder_name,'/final_save.mat'])

figure('position',[100 100 325 300],'PaperPositionMode','auto','color','white')
set(gcf,'DefaultLineLineWidth',linewidth)

axes('position',[0.23 0.15 0.72 0.8])

statistics_counter = stats.counter;
time_range = linspace(0,total_time,stats.counter);

eval([plot_name,'_plot(linewidth,stats,statistics_counter,time_range);'])

% cell_area_plot(linewidth,stats,statistics_counter,time_range);

addpath('~/Documents/export_fig/')
if save_plot_logical
   if ~exist(['Figs/',folder_name],'dir')
      mkdir('Figs/',folder_name);
   end
%    saveas(gcf,['Figs/',folder_name,'/',folder_name,'_',plot_name,'.eps'],'psc2')
   export_fig(['Figs/',folder_name,'/',folder_name,'_',plot_name,'.eps'],'-nocrop');
end

end

function mitosis_locations_plot(linewidth,stats,statistics_counter,time_range)

mitosis_locations = stats.mitosis_locations(stats.mitosis_locations(:,1)~=0,:);
mitosis_radii = sqrt(mitosis_locations(:,1).^2+mitosis_locations(:,2).^2);
if length(mitosis_radii)>1
    no_bars = 10;
    [hist_data,bin_centres] =...
        hist(mitosis_radii,no_bars);
    hist_data = hist_data./sum(hist_data)*100;
    bar(bin_centres,hist_data);
    axis([0 1.5 0 20])
end
xlabel('Mitosis radius')
ylabel('Frequency (%)')

end

function mitosis_locations_corrected_plot(linewidth,stats,statistics_counter,time_range)

mitosis_locations = stats.mitosis_locations(stats.mitosis_locations(:,1)~=0,:);
mitosis_radii = sqrt(mitosis_locations(:,1).^2+mitosis_locations(:,2).^2);
mitosis_within_original_radius = mitosis_radii(mitosis_radii<0.5);
if length(mitosis_within_original_radius)>1
    no_bars = 10;
    [hist_data,bin_centres] =...
        hist(mitosis_within_original_radius,no_bars);
    % we divide by the radius at the centre of the bin. the radius is proportional to
    % the circumference and therefore the number of cells that can fit at that radius
    hist_data = hist_data./bin_centres;
    hist_data = hist_data./sum(hist_data)*100;
    bar(bin_centres,hist_data);
    axis([0 0.5 0 20])
end
xlabel('Mitosis radius')
ylabel('Corrected frequency (%)')

end

function cell_volumes_plot(linewidth,stats,statistics_counter,time_range)

cell_volume = stats.cell_volume(1:statistics_counter,:);
plot(time_range,cell_volume(:,1))
hold on
% plot(time_range,cell_volume(:,2),'r')
% plot(time_range,cell_volume(:,3),'r')
plot(time_range,cell_volume(:,1)+2*cell_volume(:,4),'r')
plot(time_range,cell_volume(:,1)-2*cell_volume(:,4),'r')
% plot(time_range,target_volume*ones(size(time_range)),'g--')
set(gca,'FontName','arial','fontweight','bold','fontsize',13);
xlabel('Time')
ylabel('Cell volume')

axis([0 100 0 0.003])

% set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))

set(gca,'YTick',[0 0.001 0.002 0.003])
set(gca,'YTickLabel',sprintf('%0.3f|',str2num(get(gca,'YTickLabel'))*0.001))

% axes('position',[0.47 0.25 0.4 0.4])
% plot(time_range,cell_volume(:,1))
% hold on
% % plot(time_range,cell_volume(:,2),'r')
% % plot(time_range,cell_volume(:,3),'r')
% plot(time_range,cell_volume(:,1)+2*cell_volume(:,4),'r')
% plot(time_range,cell_volume(:,1)-2*cell_volume(:,4),'r')
% set(gca,'FontName','arial','fontweight','bold','fontsize',12);
% 
% axis([0 10 0 0.003])
% 
% set(gca,'YTick',[])

end

function height_to_area_ratios_plot(linewidth,stats,statistics_counter,time_range)

cell_height_to_area = stats.cell_height_to_area(1:statistics_counter,:);
plot(time_range,cell_height_to_area(:,1))
% hold on
% plot(time_range,cell_height_to_area(:,1)+2*cell_height_to_area(:,4),'r')
% plot(time_range,cell_height_to_area(:,1)-2*cell_height_to_area(:,4),'r')

set(gca,'FontName','arial','fontweight','bold','fontsize',13);
xlabel('Time')
ylabel('Height-to-area ratio')

axis([0 100 0 40])
% 
% set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
% set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))

end

function angle_deviations_plot(linewidth,stats,statistics_counter,time_range)

angle_deviation = stats.angle_deviation(1:statistics_counter,:);
plot(time_range,angle_deviation(:,1))
set(gca,'FontName','arial','fontweight','bold','fontsize',13);
xlabel('Time')
ylabel('Angle deviation')

axis([0 1 0 0.4])

set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
% set(gca,'YTickLabel',sprintf('%0.1f|',str2num(get(gca,'YTickLabel'))))

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

function cell_perimeter_plot(linewidth,stats,statistics_counter,time_range)

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

function boundary_forces_plot(linewidth,stats,statistics_counter,time_range)

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
legend boxoff

axis([0 1 0 0.1])

set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))

end

function forces_plot(linewidth,stats,statistics_counter,time_range)

total_area_force = stats.total_area_force(1:statistics_counter,:);
total_deformation_force = stats.total_deformation_force(1:statistics_counter,:);
total_elongation_force = stats.total_elongation_force(1:statistics_counter,:);
total_perimeter_force = stats.total_perimeter_force(1:statistics_counter,:);
total_tension_force = stats.total_tension_force(1:statistics_counter,:);
total_no_vertices = stats.total_no_vertices(1:statistics_counter,:);
a(1) = plot(time_range,total_area_force./total_no_vertices);
hold all
a(2) = plot(time_range,total_deformation_force./total_no_vertices);
a(3) = plot(time_range,total_elongation_force./total_no_vertices);
a(4) = plot(time_range,total_perimeter_force./total_no_vertices);
a(5) = plot(time_range,total_tension_force./total_no_vertices);
set(gca,'FontName','arial','fontweight','bold','fontsize',13);
xlabel('Time')
ylabel('Mean force per vertex')

axis([0 1 0 0.03])

set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
set(gca,'YTickLabel',sprintf('%0.3f|',str2num(get(gca,'YTickLabel'))))

legend(gca,a(1:3),'C_A','C_D','C_H','location','north')
legend boxoff
ah=axes('position',get(gca,'position'),'visible','off');
legend(ah,a(4:5),'C_P','C_L','location','northeast')
legend boxoff
set(gca,'FontName','arial','fontweight','bold','fontsize',13);

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

axis([0 100 0.002 0.016])

% set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
set(gca,'YTickLabel',sprintf('%0.3f|',str2num(get(gca,'YTickLabel'))*0.001))

end



