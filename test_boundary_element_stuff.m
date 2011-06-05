clear all; close all;

cell_dynamics
create_boundary_element

boundary_constant = 100;

figure
for i = 1:1000
	apply_boundary_force
% 	clf
% 	patchAS(node_positions(boundary_element,:))
% 	pause(0.1)
end

patchAS(node_positions(boundary_element,:))