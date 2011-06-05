disp('busy');close all;clear all;

addpath('../../export_fig/');

load Saves/generic_chemical_dependent_growth_2/final_save.mat
statistics_counter = stats.counter;
time_range = linspace(delta_t,time,statistics_counter);

figure('position',[100 100 300 300],'PaperPositionMode','auto','color','white')
axes('position',[0.18 0.18 0.75 0.75])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mitosis_locations = stats.mitosis_locations(stats.mitosis_locations(:,1)~=0,:);
% mitosis_radii = sqrt(mitosis_locations(:,1).^2+mitosis_locations(:,2).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% no_bars = 10;
% [hist_data,bin_centres] =...
%     hist(mitosis_radii,no_bars);
% hist_data = hist_data/length(mitosis_locations)*100;
% bar(bin_centres,hist_data);
% xlabel('Radius')
% ylabel('Frequency (%)')
% axis([0 1.0 0 20])
% set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))

% no_bars = 10;
% [hist_data,bin_centres] = hist(mitosis_locations(:,1),no_bars);
% hist_data = hist_data/length(mitosis_locations)*100;
% bar(bin_centres,hist_data);
% xlabel('x-Location')
% ylabel('Frequency (%)')
% set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
% axis([-1.0 1.0 0 25])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('outerposition',[100 100 300 300],'PaperPositionMode','auto')
% mitosis_radii_within_original_radius = mitosis_radii(mitosis_radii<0.5);
% no_bars = 10;
% [hist_data,bin_centres] =...
% 	hist(mitosis_radii_within_original_radius,no_bars);
% hist_data = hist_data./bin_centres;
% hist_data = hist_data./sum(hist_data)*100;
% hist_data = hist_data/length(mitosis_locations)*100;
% radii_bin_centres = linspace(0,0.5,11);
% bar(bin_centres,hist_data);
% xlabel('Radius')
% ylabel('Frequency (%)')
% axis([0 0.5 0 20])
% set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))

% figure('outerposition',[100 100 300 300],'PaperPositionMode','auto')
% mitosis_within_original_radius = mitosis_locations(mitosis_radii<0.5,1);
% no_bars = 10;
% [hist_data,bin_centres] =...
% 	hist(mitosis_within_original_radius,no_bars);
% hist_data = hist_data./(2*sqrt(0.5^2 - bin_centres.^2));
% hist_data = hist_data./sum(hist_data)*100;
% bar(bin_centres,hist_data);
% xlabel('x-Location')
% ylabel('Frequency (%)')
% axis([-0.5 0.5 0 25])
% set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% total_Dpp = stats.total_Dpp(1:statistics_counter,:);
% plot(time_range,total_Dpp,'linewidth',3)
% axis([0 25000 0.05 0.09])
% 
% % set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
% set(gca,'XTickLabel',['0.0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'])
% set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))
% 
% xlabel('Time (normalised)')
% ylabel('Total morphogen quantity')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% triangle_quality_stats = stats.triangle_quality(1:statistics_counter,:);
% plot(time_range,triangle_quality_stats(:,1))
% hold on
% plot(time_range,triangle_quality_stats(:,1)+2*triangle_quality_stats(:,4),'r')
% axis([0 25000 0 20])
% set(gca,'XTickLabel',['0.0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'])
% 
% xlabel('Time (normalised)')
% ylabel('Triangle quality')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xy_ratio = stats.xy_ratio(1:statistics_counter,:);
% plot(time_range,xy_ratio)
% xlabel('Time (normalised)')
% ylabel('XY ratio')
% axis([0 25000 0.9 1.1])
% set(gca,'XTickLabel',['0.0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'])
% set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gradient_volume_distribution_x = stats.normalised_volume_distribution_x(:,1);
% plot(time_range,gradient_volume_distribution_x)
% xlabel('Time')
% ylabel('Volume distribution gradient (x)')
% axis([0 1 -0.2 0.2])
% set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
% set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))

% gradient_volume_distribution_y = stats.normalised_volume_distribution_y(:,1);
% plot(time_range,gradient_volume_distribution_y)
% xlabel('Time')
% ylabel('Volume distribution gradient (y)')
% axis([0 1 -0.2 0.2])
% set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
% set(gca,'YTickLabel',sprintf('%0.1f|',str2num(get(gca,'YTickLabel'))))

% gradient_area_distribution_x = stats.normalised_area_distribution_x(:,1);
% plot(time_range,gradient_area_distribution_x)
% xlabel('Time')
% ylabel('Area distribution gradient (x)')
% axis([0 1 -0.2 0.2])
% set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
% set(gca,'YTickLabel',sprintf('%0.1f|',str2num(get(gca,'YTickLabel'))))

% export_fig Figs/generic_chemical_dependent_growth_2/volume_distribution_x.jpg -r300 -nocrop

% mean_chemical_ingestion_rate = stats.chemical_ingestion_rate(1:statistics_counter,1);
% % max_chemical_internalised = stats.chemical_internalised(1:statistics_counter,2);
% % min_chemical_internalised = stats.chemical_internalised(1:statistics_counter,3);
% plot(time_range,mean_chemical_ingestion_rate)
% axis([0 1 0 12e-7])
% hold on
% % plot(time_range,max_chemical_internalised,'r')
% % plot(time_range,min_chemical_internalised,'r')
% xlabel('Time')
% ylabel('Mean chemical ingestion rate')
% set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
% % set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))
% 
% export_fig Figs/generic_chemical_dependent_growth_2/mean_chemical_ingested.jpg -r300 -nocrop

total_chemical_released = stats.chemical_source(1:statistics_counter,1);
total_internal_chemical = stats.total_internal_chemical(1:statistics_counter,1);
net_chemical_released = total_chemical_released - total_internal_chemical;
% subplot(2,3,2)
plot(time_range,total_chemical_released)
hold all
plot(time_range,total_internal_chemical)
plot(time_range,net_chemical_released)
xlabel('Time')
ylabel('Net chemical released')
legend('Source','Ingested','Net','Location','Northwest')
if max(total_chemical_released)>0
    axis([0 time 0.9*min(total_chemical_released) 1.1*max(total_chemical_released)])
else
    axis_values = axis;
    axis([axis_values(1) time axis_values(3) axis_values(4)])
end

