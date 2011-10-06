disp('busy'); clear all; close all;

saved_iterations = [1 200:200:1000];

for i = 1:length(saved_iterations)
   
   current_saved_iteration = saved_iterations(i);
   load(['Saves/drosophila_epidermis/iteration_',num2str(current_saved_iteration)]);
   
   cell_mean_x_positions = cellfun(@(x)mean(vertices.position(x,1)),cells.vertices);

   shifted_cell_mean_x_positions = cell_mean_x_positions-min(cell_mean_x_positions);
   normalised_cell_mean_x_positions = shifted_cell_mean_x_positions/max(shifted_cell_mean_x_positions)-0.5;
   
   cells_logical = cells.state~=3 & cells.state~=4;
   
%    figure;
%    plot(normalised_cell_mean_x_positions(cells_logical),cells.volume(cells_logical),'o');
% 
%    axis([-0.5 0.5 0 1e-3])
   
   figure;
   plot(normalised_cell_mean_x_positions(cells_logical),cells.internal_chemical_quantity(cells_logical),'o')
   
   axis([-0.5 0.5 0 0.35])
   
end
%    a = polyfit(cell_centre_positions_x,cells.volume,1);
% hold on;
% plot(cell_centre_positions_x,a(1)*cell_centre_positions_x+a(2))