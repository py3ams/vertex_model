disp('busy');clear all; close all;

refinement_0 = load('Saves/refinement_0');
refinement_1 = load('Saves/refinement_1');
refinement_2 = load('Saves/refinement_2');
refinement_3 = load('Saves/refinement_3');

length_vertices = length(refinement_0.vertices.position);
length_cells = length(refinement_0.cells.vertices);

norms(1) = sqrt(sum([refinement_0.FEM_nodes.concentration([1:length_vertices ...
    end-length_cells:end])-refinement_1.FEM_nodes.concentration([1:length_vertices ...
    end-length_cells:end])].^2));
    
norms(2) = sqrt(sum([refinement_1.FEM_nodes.concentration([1:length_vertices ...
    end-length_cells:end])-refinement_2.FEM_nodes.concentration([1:length_vertices ...
    end-length_cells:end])].^2));

norms(3) = sqrt(sum([refinement_2.FEM_nodes.concentration([1:length_vertices ...
    end-length_cells:end])-refinement_3.FEM_nodes.concentration([1:length_vertices ...
    end-length_cells:end])].^2));

plot(norms)