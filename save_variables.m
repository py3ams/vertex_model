function save_variables(filename,cells,vertices,FEM_elements,FEM_nodes,cell_growth_speeds_matrix,refined_edge_matrix)

eval(['save ',filename,' cells vertices FEM_elements FEM_nodes cell_growth_speeds_matrix refined_edge_matrix']);