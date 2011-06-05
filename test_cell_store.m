function broken_nodes = test_cell_store(cell_store,cells)

broken_nodes = [];

cells_logical = ~cellfun('isempty',cells);

[dummy_cell_store] =...
    CreateCellStore(cells(cells_logical),length(cell_store));
   
logical_equal = sum(cell_store,2)==sum(dummy_cell_store,2);

if ~all(logical_equal)
    broken_nodes = find(~logical_equal);
    error(['**** cell_store broken in ',function_name,' at iteration ',...
        int2str(iteration),' ****'])
end

