% can use this file to draw all the statistical plots from a simulation if
% for some reason you have the saves but not the figs

clear all; close all; disp('busy');

simulation_name = 'morphogen_dependent_growth';

% load(['Saves/',simulation_name,'/final_save']);
load(['Saves/',simulation_name,'/iteration_30000']);

fig_saves_logical = true;

screen_size = get(0,'screensize');
figure_position =...
    [(screen_size(3)-5/4*screen_size(4))/2 0 5/4*screen_size(4) screen_size(4)];

% visualiser(cells,node_positions,Brk,Dpp,FEM_elements,previous_FEM_node_positions,...
%     axis_values,axis_values_FEM,FEM_angle,figure_position,false,...
%     0,movie_logical,shading_style,update_period_1,update_period_2,...
%     view_FEM_logical,visualise_initial_configuration);
 
 visualiser(cells,node_positions,Brk,Dpp,FEM_elements,previous_FEM_node_positions,...
    axis_values,axis_values_FEM,FEM_angle,figure_position,false,...
    0,movie_logical,no_original_cells,original_cells_colour,shading_style,...
    update_period_1,update_period_2,view_FEM_logical,visualise_initial_configuration);

% visualiser(cells,node_positions,Brk,Dpp,FEM_elements,previous_FEM_node_positions,...
%     axis_values,axis_values_FEM,FEM_angle,figure_position,true,...
%     iteration,movie_logical,shading_style,update_period_1,update_period_2,...
%     view_FEM_logical,visualise_initial_configuration);

visualiser(cells,node_positions,Brk,Dpp,FEM_elements,previous_FEM_node_positions,...
    axis_values,axis_values_FEM,FEM_angle,figure_position,false,...
    iteration,movie_logical,no_original_cells,original_cells_colour,shading_style,...
    update_period_1,update_period_2,view_FEM_logical,visualise_initial_configuration);

statistical_plots(delta_t,fig_saves_location,fig_saves_logical,figure_position,...
    movie_location,movie_logical,movie_name,no_frames_for_statistical_plots,...
    statistics_counter,stats,time);
