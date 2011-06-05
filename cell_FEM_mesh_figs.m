disp('busy');tic;close all; clear all;%profile on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% s = RandStream.create('mt19937ar','seed',5489);
% s = RandStream.create('mt19937ar','seed',sum(100*clock));
% RandStream.setDefaultStream(s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_iterations = 500;
simulation_name = '';

target_no_cells_logical = false;
target_no_cells = 400;

delta_t = 1;
viscosity = 1;

fig_saves_logical = false;
fig_saves_name = simulation_name;

full_saves_logical = false;
full_saves_name = simulation_name;
full_saves_period = max(floor(max_iterations/5),1);

regular_tests_logical = false;

statistics_period = max(floor(max_iterations/1000),1);

%%%%%%%%%%%%%%%%%%%%%%%%% Initial configuration parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

anneal_initial_configuration_logical = true;
array_sizes = 8000;
compile_mex_functions_logical = false;
configuration_noise = 0.8;
file_to_load = 'initial_save';
grid_size = [7,7];
load_from_file_logical = false;
no_original_cells = prod(grid_size);
original_cells_colour = 'r';

if target_no_cells_logical && prod(grid_size) > target_no_cells
    error('Target number of cells is smaller than initial number');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update position parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

update_positions_logical = true;
update_positions_start = 0;

cell_volume_growth_logical = true;
cell_volume_growth_type = 'logistic';

average_cell_growth_speeds = [0.01 0.02];
% no_growth_time = 5000;
no_growth_time = 0;
% average_cell_growth_speeds = [5e-7 1e-6];
% medial_lateral_threshold_factor = 0.5;
medial_lateral_threshold_factor = 100;
% medial_lateral_threshold_factor = 0.0;
target_area_factor = 1.0;
target_volume_factor = 1.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Vertex rearrangement parameters %%%%%%%%%%%%%%%%%%%%%%%%%

T1_swaps_logical = true;
T1_swaps_start = 0;
T1_probability = 1.0;

threshold_T1_swaps_factor = 0.02;

% protection_time = 100;
protection_time = 0;

mitosis_logical = true;
mitosis_start = 0;

mitosis_angles_type = 'uniform';
% mitosis_angles = [0 0];
mitosis_dependence = 'volume';
mitosis_random_logical = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEM parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FEM_solve_logical = true;

% degradation_constant = 0.005;
degradation_constant = 0;
diffusion_speed = 0.001;
% source_magnitude = 0.05;
source_magnitude = 0;
source_width = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Force constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

area_force_constant = 1e-0;
deformation_force_constant = 1e-4;
elongation_force_constant = 5e-7;
perimeter_force_constant = 1e-4;
tension_force_constant = 2e-4;

% area_force_constant = 0;
% deformation_force_constant = 0;
% elongation_force_constant = 0;
% perimeter_force_constant = 0;
% tension_force_constant = 0;

boundary_deformation_force_constant = 1e-3;
boundary_edge_force_constant = 1e-1;

% boundary_deformation_force_constant = 0;
% boundary_edge_force_constant = 0;

tension_anisotropy_factor = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Movie parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

movie_logical = 0;

% axis_values = [-1 1 -1 1];
axis_values = 'equal';
% axis_values_FEM = [-1 1 -1 1 -0.5 1.5];
% axis_values_FEM = [-1 1 -1 1 -0.1 0.2];
axis_values_FEM = 'equal';
extra_pause = 0.0;
% extra_pause = 0.1;
FEM_angle = [45 45];
movie_name = simulation_name;
movie_period = 50;
movie_start = 0;
no_frames_for_statistical_plots = 100;
shading_style = 'faceted';
update_period_1 = 1e9;
update_period_2 = movie_period;
view_FEM_logical = true;
visualise_initial_configuration = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compile mex functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compile_mex_functions(compile_mex_functions_logical);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate initial configuration %%%%%%%%%%%%%%%%%%%%%%%%%%

[cells,node_positions,cell_elements,cell_growth_speeds,...
    cell_growth_speeds_matrix,cell_store,cell_volumes,cells_per_node,Dpp,...
    FEM_elements,FEM_node_positions,mitosis_angles,previous_FEM_node_positions] =...
    initial_configuration(anneal_initial_configuration_logical,array_sizes,...
    average_cell_growth_speeds,configuration_noise,grid_size,FEM_solve_logical,...
    file_to_load,load_from_file_logical,medial_lateral_threshold_factor,...
    mitosis_angles_type,target_no_cells,area_force_constant,...
    boundary_deformation_force_constant,boundary_edge_force_constant,...
    deformation_force_constant,elongation_force_constant,perimeter_force_constant,...
    tension_force_constant);


figure('PaperPositionMode','auto','outerposition',[100 100 900 900])
axes('position',[0 0 1 1])
figure_loop(cells,node_positions,length(cells),[255,30,30]/255,5)
axis equal
axis off

figure('PaperPositionMode','auto','outerposition',[100 100 900 900])
axes('position',[0 0 1 1])
cellfun(@(x)patch(node_positions(x,1),node_positions(x,2),[255,30,30]/255,'linewidth',5,'FaceAlpha',0),cells);
axis equal	
axis off
grid off
hold on
trisurf(FEM_elements,FEM_node_positions(:,1),FEM_node_positions(:,2),zeros(length(FEM_node_positions),1),'linewidth',1.5)
view(0,90)
colormap([50,180,50]/255)
