disp('busy');clear all; close all;

mex T1Swaps.cpp

before = load('T1_swaps_save_2');

before.threshold_T1_swaps = 0.5*before.mean_edge_length;

after.cells = before.cells;
after.vertices = before.vertices;

[after.cells.vertices,after.vertices.position,after.cells.FEM_elements,after.vertices.cells,after.vertices.no_cells,after.FEM_nodes.concentration,...
    after.FEM_elements.nodes,no_T1_swaps_this_iteration,after.FEM_nodes.previous_position,...
    after.vertices.time_created] = T1Swaps(before.cells.vertices,before.vertices.position,before.cells.FEM_elements,...
    before.vertices.cells,before.vertices.no_cells,before.FEM_nodes.concentration,before.FEM_elements.nodes,before.FEM_solve_logical,...
    before.FEM_nodes.previous_position,before.protection_time,before.T1_probability,...
    before.threshold_T1_swaps,before.time,before.vertices.time_created);

figure('units','normalized','outerposition',[0 0 1 0.5])

subplot(1,2,1)
figure_loop(before.cells,before.vertices)
% axis([-0.2 0.1 -0.1 0.1])

subplot(1,2,2)
figure_loop(after.cells,after.vertices)
% axis([-0.2 0.1 -0.1 0.1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('units','normalized','outerposition',[0 0 1 1])
FEM_element_nodes = before.FEM_elements.nodes(before.FEM_elements.nodes(:,1)>0,:);

subplot(2,2,1)
trisurf(FEM_element_nodes,before.FEM_nodes.previous_position(:,1),...
    before.FEM_nodes.previous_position(:,2),zeros(length(before.FEM_nodes.previous_position),1))
hold on;
cellfun(@(x)patch(before.vertices.position(x,1),before.vertices.position(x,2),'r','linewidth',2.5,'FaceAlpha',0),before.cells.vertices);
% axis([-0.2 0.1 -0.1 0.1 0 0.4])
view(0,90)

FEM_element_nodes = after.FEM_elements.nodes(after.FEM_elements.nodes(:,1)>0,:);

subplot(2,2,2)
trisurf(FEM_element_nodes,after.FEM_nodes.previous_position(:,1),...
    after.FEM_nodes.previous_position(:,2),zeros(length(after.FEM_nodes.previous_position),1))
hold on;
cellfun(@(x)patch(after.vertices.position(x,1),after.vertices.position(x,2),'r','linewidth',2.5,'FaceAlpha',0),after.cells.vertices);
% axis([-0.2 0.1 -0.1 0.1 0 0.4])
view(0,90)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('units','normalized','outerposition',[0 0 1 1])
FEM_element_nodes = before.FEM_elements.nodes(before.FEM_elements.nodes(:,1)>0,:);

subplot(2,2,3)
trisurf(FEM_element_nodes,before.FEM_nodes.previous_position(:,1),...
    before.FEM_nodes.previous_position(:,2),before.FEM_nodes.concentration)
% axis([-0.2 0.1 -0.1 0.1 0 0.4])
view(-7,36)

FEM_element_nodes = after.FEM_elements.nodes(after.FEM_elements.nodes(:,1)>0,:);

subplot(2,2,4)
trisurf(FEM_element_nodes,after.FEM_nodes.previous_position(:,1),...
    after.FEM_nodes.previous_position(:,2),after.FEM_nodes.concentration)
% axis([-0.2 0.1 -0.1 0.1 0 0.4])
view(-7,36)