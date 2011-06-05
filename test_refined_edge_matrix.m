% disp('busy')

no_FEM_nodes = length(FEM_nodes.edge);

no_refined_edges = 0;

for current_FEM_node = 1:no_FEM_nodes
	
	if FEM_nodes.edge(current_FEM_node)>0
		
		no_refined_edges = no_refined_edges+1;
		
		if ~refined_edge_matrix(FEM_nodes.edge(current_FEM_node,1),FEM_nodes.edge(current_FEM_node,2)) ||...
				~refined_edge_matrix(FEM_nodes.edge(current_FEM_node,2),FEM_nodes.edge(current_FEM_node,1))
			
			error(['warning : refined_edge_matrix missing an entry for FEM_node ',num2str(current_FEM_node)]);
			
		end
		
	end
	
end

if sum(sum(refined_edge_matrix))>2*no_refined_edges
	
	error('warning : too many entries in refined_edge_matrix')
	
end

