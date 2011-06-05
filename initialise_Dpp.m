function Dpp = initialise_Dpp(cells_per_node,FEM_node_positions,no_cells)

no_FEM_nodes = length(FEM_node_positions);

Dpp = zeros(no_FEM_nodes,1);

width = 0.15;

for i = 1:length(cells_per_node)
    
%     if cells_per_node(i) > 0 && abs(FEM_node_positions(i,1)) < width
        
    if cells_per_node(i) > 0 && sqrt((FEM_node_positions(i,1).^2)+(FEM_node_positions(i,2).^2))<width
		 
        Dpp(i) = 1;
        
    end
    
%     if cells_per_node(i) > 0
%         
%         Dpp(i) = exp(-abs(10*FEM_node_positions(i,1)))*(1 + 0.2*randn);
%         
%     end
    
end

for i = (no_FEM_nodes-no_cells):no_FEM_nodes
    
    if sqrt((FEM_node_positions(i,1).^2)+(FEM_node_positions(i,2).^2))<width
        
        Dpp(i) = 1;
        
    end
    
%     Dpp(i) = exp(-abs(10*FEM_node_positions(i,1)))*(1 + 0.2*randn);
    
end