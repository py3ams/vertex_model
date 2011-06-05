function [Brk,Dpp] = initialise_morphogens(cells_per_node,FEM_node_positions)

no_FEM_nodes = length(FEM_node_positions);
no_regular_nodes = length(cells_per_node);
no_cells = no_FEM_nodes - no_regular_nodes;

Brk = zeros(no_FEM_nodes,1);
Dpp = zeros(no_FEM_nodes,1);

for i = 1:length(cells_per_node)
    
    if cells_per_node(i) > 0
        
        if abs(FEM_node_positions(i,1)) < 0.0
            
            Dpp(i) = 1;
            
        else
            Brk(i) = 1;
            
        end
        
        %     if cells_per_node(i) > 0
        %
        %         Dpp(i) = exp(-abs(10*FEM_node_positions(i,1)))*(1 + 0.2*randn);
        %
        %     end
        
    end
    
end

for i = 1:no_cells
    
    if abs(FEM_node_positions(i+no_regular_nodes,1)) < 0.0
        
        Dpp(i+no_regular_nodes) = 1;
        
    else
        
        Brk(i+no_regular_nodes) = 1;
        
    end
    
    %     Dpp(i) = exp(-abs(10*FEM_node_positionse_positions(i,1)))*(1 + 0.2*randn);
    
end