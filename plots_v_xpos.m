disp('busy'); clear all; close all;

save_plot_logical = true;
initial_plot_logical = true;

saved_iterations = [1000:1000:5000];

folder_name = 'drosophila_epidermis_limited_spi';
plot_name = 'internal_chemical_quantity';

alphabet = 'abcdefghijklmnopqrstuvwxyz';
alphabet_index = 1;

linewidth = 2;

if initial_plot_logical
   
   load(['Saves/',folder_name,'/initial_save']);

   cell_mean_x_positions = cellfun(@(x)mean(vertices.position(x,1)),cells.vertices);

   shifted_cell_mean_x_positions = cell_mean_x_positions-min(cell_mean_x_positions);
   normalised_cell_mean_x_positions = shifted_cell_mean_x_positions/max(shifted_cell_mean_x_positions)-0.5;
   
   cells_logical = cells.state~=3 & cells.state~=4;
      
   figure('position',[100 100 325 300],'PaperPositionMode','auto','color','white')
   set(gcf,'DefaultLineLineWidth',linewidth)
   axes('position',[0.23 0.15 0.72 0.8])
   
   plot_style = 'o';
   eval(['plot(normalised_cell_mean_x_positions(cells_logical),cells.',plot_name,'(cells_logical),plot_style)'])
   xlim([-0.5 0.5])
   ylim([0 0.25])
   set(gca,'FontName','arial','fontweight','bold','fontsize',13);
   
   set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
   set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))
   
   addpath('~/Documents/export_fig/')
   if save_plot_logical
      if ~exist(['Figs/',folder_name],'dir')
         mkdir('Figs/',folder_name);
      end
      %    saveas(gcf,['Figs/',folder_name,'/',folder_name,'_',plot_name,'.eps'],'psc2')
      export_fig(['Figs/',folder_name,'/',folder_name,'_',plot_name,'_',alphabet(alphabet_index),'.eps'],'-nocrop');
      alphabet_index = alphabet_index+1;
   end
   
end
   
for i = 1:length(saved_iterations)
   
   current_saved_iteration = saved_iterations(i);
   load(['Saves/',folder_name,'/iteration_',num2str(current_saved_iteration)]);
   
   cell_mean_x_positions = cellfun(@(x)mean(vertices.position(x,1)),cells.vertices);

   shifted_cell_mean_x_positions = cell_mean_x_positions-min(cell_mean_x_positions);
   normalised_cell_mean_x_positions = shifted_cell_mean_x_positions/max(shifted_cell_mean_x_positions)-0.5;
   
   cells_logical = cells.state~=3 & cells.state~=4;
      
   figure('position',[100 100 325 300],'PaperPositionMode','auto','color','white')
   set(gcf,'DefaultLineLineWidth',linewidth)
   axes('position',[0.23 0.15 0.72 0.8])
   
   plot_style = 'o';
   eval(['plot(normalised_cell_mean_x_positions(cells_logical),cells.',plot_name,'(cells_logical),plot_style)'])
   xlim([-0.5 0.5])
   ylim([0 0.25])
   set(gca,'FontName','arial','fontweight','bold','fontsize',13);
   
   xlabel('x-position')
   ylabel('Internal chemical')
   
   set(gca,'XTickLabel',sprintf('%0.1f|',str2num(get(gca,'XTickLabel'))))
   set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))

   if save_plot_logical
      export_fig(['Figs/',folder_name,'/',folder_name,'_',plot_name,'_',alphabet(alphabet_index),'.eps'],'-nocrop');
      alphabet_index = alphabet_index+1;
   end
   
end
