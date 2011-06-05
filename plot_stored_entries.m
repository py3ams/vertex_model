disp('busy');%close all

no_stored_entries  = 3;

for i = 1:no_stored_entries
    
    figure
    
    for j =1:length(cells.previous_vertices{i})
        hold on
        patchAS(vertices.previous_positions(cells.previous_vertices{i}{j},2*i-1:2*i),'r',1)
    end
    
    
%     figure_loop(cells.previous_vertices{i},vertices.previous_position(:,2*i-1:2*i))
%     axis([-0.5403 -0.4274 -0.2500 -0.1915])
%     axis tight
    
end